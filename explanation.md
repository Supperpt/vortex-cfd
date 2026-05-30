# Como funciona o vortex-cfd

Este documento explica o que o programa faz, como o OpenFOAM funciona por baixo, o que significam os outputs actuais (U e p), e como obter outros outputs relevantes para o doutoramento (WSS, gradiente de WSS, velocidade num ponto, etc.).

---

## 1. Visão geral do pipeline

```
STLs do VMTK
    │
    ▼
[vortex-cfd: Python]
    │  · Detecta OpenFOAM
    │  · Pede ao utilizador para identificar wall / inlet / outlet
    │  · Escala mm → m
    │  · Calcula parâmetros geométricos (bounding box, área do inlet, ponto interior)
    │  · Gera todos os ficheiros de configuração do OpenFOAM (13 ficheiros)
    │
    ▼
[surfaceFeatureExtract]   Extrai as arestas das superfícies STL (ficheiros .eMesh)
    │
    ▼
[blockMesh]               Cria uma malha hexaédrica simples em volta do vaso (o "cubo")
    │
    ▼
[snappyHexMesh]           Refina e encaixa a malha na geometria do vaso
    │  · Fase 1 (castellated): corta células que intersectam a parede do vaso
    │  · Fase 2 (snap): encaixa os vértices à superfície STL
    │  · Fase 3 (addLayers): adiciona 4 camadas prismáticas junto à parede
    │
    ▼
[checkMesh]               Verifica qualidade da malha (ortogonalidade < 70°, skewness < 4)
    │
    ▼
[decomposePar]            Divide a malha em N partes para correr em paralelo
    │
    ▼
[pimpleFoam]              Resolve as equações de Navier-Stokes (o solver CFD)
    │  · Lê as condições iniciais e de fronteira
    │  · Avança no tempo passo a passo (ciclo cardíaco pulsátil)
    │  · Escreve snapshots a cada T/50 (50 por ciclo)
    │
    ▼
[reconstructPar]          Junta os resultados dos N processos num só directório
    │
    ▼
Directório case_YYYYMMDD_HHMMSS/
    └── Abre em ParaView
```

---

## 2. O que é a malha e porque é importante

O OpenFOAM não trabalha com a geometria STL directamente — precisa de uma **malha volumétrica**: o interior do vaso dividido em milhares de pequenas células (hexaedros e prismas). As equações de Navier-Stokes são resolvidas em cada célula.

### blockMesh — a malha de fundo

O `blockMesh` cria um cubo simples em volta do vaso com células de ~2 mm. Este cubo é o ponto de partida — ainda não tem a forma do vaso.

### snappyHexMesh — a malha final

O `snappyHexMesh` pega no cubo e:
1. Remove as células fora do vaso (o espaço exterior)
2. Refina as células perto da parede (até ~0.125 mm)
3. Encaixa os vértices na superfície STL
4. Adiciona 4 camadas prismáticas junto à parede

As **camadas prismáticas** são críticas para o WSS. A tensão de corte na parede é calculada a partir do gradiente de velocidade dU/dr — quanto mais fina for a primeira célula junto à parede, mais preciso é esse gradiente. Sem camadas, o WSS seria uma aproximação grosseira.

### locationInMesh — o ponto que define "dentro"

O snappyHexMesh precisa de saber qual é o lado "fluido". Fornecemos um ponto que está garantidamente dentro do lúmen: o centroide do inlet deslocado um raio para dentro ao longo da normal da superfície. Se este ponto estiver errado, a malha fica invertida (simula o exterior do vaso) — foi exactamente o que aconteceu na primeira corrida.

---

## 3. O solver: pimpleFoam

O `pimpleFoam` resolve as **equações de Navier-Stokes incompressíveis** em regime transiente:

```
∂U/∂t + (U·∇)U = -∇p + ν∇²U      (equação do momento)
∇·U = 0                             (incompressibilidade)
```

Onde:
- **U** = campo de velocidade (vector 3D em cada ponto)
- **p** = pressão cinemática (p/ρ, em m²/s²)
- **ν** = viscosidade cinemática = 3.3×10⁻⁶ m²/s

### Como avança no tempo

O timestep é **adaptativo**: começa em 10⁻⁵ s e ajusta-se automaticamente para manter o número de Courant Co < 0.8. O número de Courant é Co = U·Δt/Δx — mede quantas células o fluido atravessa por timestep. Se Co > 1 o solver diverge. Na sístole a velocidade é ~10× maior que na diástole, por isso o timestep reduz automaticamente no pico sistólico.

### Condições de fronteira

| Superfície | Velocidade (U) | Pressão (p) |
|---|---|---|
| Parede (wall) | noSlip — U = 0 | zeroGradient |
| Inlet | flowRateInletVelocity — caudal pulsátil | zeroGradient |
| Outlets | inletOutlet — zero-gradient saída, bloqueia refluxo | fixedValue = 0 Pa |
| Background (cubo) | slip | zeroGradient |

O inlet usa uma tabela de (tempo, caudal) gerada a partir da forma de onda cardíaca multiplicada pela área do inlet e pela velocidade média fornecida pelo utilizador.

### Snapshots

O solver escreve 50 snapshots por ciclo cardíaco. Cada snapshot é uma pasta com o tempo (ex: `0.01714/`) que contém os ficheiros `U` e `p` com os valores em todas as células.

---

## 4. Os outputs actuais: U e p

### U — campo de velocidade

**Unidades:** m/s  
**Tipo:** vector 3D — em cada célula tens (Ux, Uy, Uz)  
**No ParaView:** podes ver a magnitude |U| ou cada componente separadamente

O que deves ver:
- **0 m/s na parede** (condição noSlip)
- **Máximo no centro** do vaso (perfil de Poiseuille)
- **Zona de recirculação** dentro do saco do aneurisma com velocidades baixas e direcções variáveis
- **Pulsação temporal** — U aumenta no pico sistólico (~30% do ciclo), diminui na diástole

Valores típicos para a artéria carótida interna: 0.1–1.0 m/s.

### p — pressão cinemática

**Unidades:** m²/s² (pressão cinética = p_real / ρ)  
**Para converter para Pa:** multiplica por ρ = 1060 kg/m³  
**Tipo:** escalar

O que deves ver:
- Gradiente do inlet para os outlets (inlet maior, outlets = 0 por definição)
- Valores típicos: 0–200 Pa (0–0.19 m²/s² em pressão cinemática)
- O outlet está fixo a 0 Pa — só importa a pressão *relativa*, não absoluta (fluido incompressível)

---

## 5. Como obter outros outputs

### 5.1 WSS — Wall Shear Stress (tensão de corte na parede)

O WSS é o output clínico mais importante. É a tensão tangencial que o sangue exerce na parede:

```
WSS = μ · (∂U/∂n)|parede
```

Onde n é a direcção normal à parede e μ = ρν = 3.71×10⁻³ Pa·s.

**Fase C do roadmap implementará isto automaticamente.** Por agora, podes calcular no ParaView:

1. Abre o caso
2. `Filters → Gradient Of Unstructured DataSet` sobre U
3. Depois `Filters → Calculator` para extrair a componente tangencial

Ou de forma mais directa, adicionar ao `controlDict` antes de correr:

```cpp
functions
{
    wallShearStress
    {
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        patches         (wall);
    }
}
```

Isto faz o pimpleFoam calcular e guardar `wallShearStress` em cada snapshot automaticamente.

**Valores típicos:** 0–50 Pa. WSS muito baixo (<0.4 Pa) e oscilatório é associado a risco de ruptura.

### 5.2 TAWSS — Time-Averaged WSS

Média do WSS ao longo de um ciclo cardíaco completo:

```
TAWSS = (1/T) ∫₀ᵀ |WSS(t)| dt
```

Adicionar ao `controlDict`:

```cpp
fieldAverage1
{
    type            fieldAverage;
    libs            (fieldFunctionObjects);
    writeControl    writeTime;
    timeStart       0.857;   // começa no 2º ciclo (descarta o 1º como transiente)
    fields
    (
        wallShearStress { mean on; prime2Mean off; base time; }
    );
}
```

Isto gera `wallShearStressMean` — a média temporal do WSS na parede.

### 5.3 OSI — Oscillatory Shear Index

Mede o quanto a direcção do WSS oscila durante o ciclo:

```
OSI = 0.5 × (1 - |∫WSS dt| / ∫|WSS| dt)
```

- OSI = 0: WSS sempre na mesma direcção (fluxo unidireccional, saudável)
- OSI = 0.5: WSS completamente oscilatório (recirculação, risco de ruptura)

Não existe function object nativo para OSI no OpenFOAM — requer pós-processamento. A Fase C vai implementar isto em Python com pyvista, lendo os snapshots do WSS.

### 5.4 Velocidade num ponto específico do aneurisma

**No ParaView:** `Filters → Probe Location` — defines as coordenadas X, Y, Z de um ponto e obtens U e p nesse ponto ao longo do tempo. Exportas para CSV com `File → Save Data`.

**Pelo OpenFOAM** (antes de correr), adicionas ao `controlDict`:

```cpp
probes
{
    type        probes;
    libs        (sampling);
    writeControl timeStep;
    writeInterval 1;
    probeLocations
    (
        (0.012 0.008 0.005)   // coordenadas em metros do ponto a amostrar
    );
    fields (U p);
}
```

Gera um ficheiro `postProcessing/probes/0/U` com a série temporal.

### 5.5 Perfil de velocidade numa secção transversal

**No ParaView:** `Filters → Slice` sobre o internalMesh, posicionas o plano no local desejado, depois `Filters → Plot Over Line` para ver o perfil radial.

**Pelo OpenFOAM**, adicionas ao `controlDict`:

```cpp
surfaces
{
    type        surfaces;
    libs        (sampling);
    writeControl writeTime;
    surfaceFormat vtk;
    fields (U p);
    surfaces
    (
        slice_aneurysm
        {
            type        cuttingPlane;
            point       (0.012 0.008 0.005);
            normal      (1 0 0);
            interpolate true;
        }
    );
}
```

### 5.6 Caudal (flow rate) numa secção

Para validar que o caudal simulado é fisiológico:

```cpp
flowRatePatch
{
    type        surfaceFieldValue;
    libs        (fieldFunctionObjects);
    writeControl writeTime;
    regionType  patch;
    name        inlet;
    operation   sum;
    fields      (phi);   // phi é o fluxo volumétrico por face
}
```

---

## 6. Como adicionar estes outputs à simulação actual

Para qualquer um dos outputs acima, editas o ficheiro `system/controlDict` dentro da pasta do caso **antes** de correr o pimpleFoam, adicionando o bloco `functions { ... }` no fim.

Se já correste e queres recalcular WSS nos resultados existentes sem re-simular:

```bash
cd /tmp/vortex_test/case_20260530_XXXXXX
# Adiciona wallShearStress ao controlDict, depois:
pimpleFoam -postProcess -func wallShearStress
```

O `-postProcess` lê os snapshots já calculados e aplica a function object sem re-simular.

---

## 7. Estrutura do directório de resultados

```
case_20260530_XXXXXX/
├── case_20260530_XXXXXX.foam   ← abre no ParaView
├── 0/                          ← condições iniciais (t=0)
│   ├── U                       ← velocidade inicial (zero)
│   └── p                       ← pressão inicial (zero)
├── 0.01714/                    ← snapshot t=0.01714 s (1º de 50 do ciclo)
│   ├── U                       ← campo de velocidade neste instante
│   └── p                       ← campo de pressão neste instante
├── ... (mais snapshots)
├── constant/
│   ├── triSurface/             ← os STLs copiados (wall.stl, inlet.stl, ...)
│   ├── transportProperties     ← ρ=1060, ν=3.3e-6
│   └── turbulenceProperties    ← laminar
├── system/
│   ├── controlDict             ← duração, timestep, writeInterval
│   ├── fvSchemes               ← esquemas numéricos (backward, linearUpwind)
│   ├── fvSolution              ← tolerâncias do solver linear
│   ├── snappyHexMeshDict       ← parâmetros de refinamento da malha
│   └── ...
└── patch_labels.json           ← mapeamento nome-STL → tipo (wall/inlet/outlet)
```

Cada pasta com um número (ex: `0.01714/`) é um snapshot — o OpenFOAM escreve 50 por ciclo cardíaco. Quando abres o `.foam` no ParaView, ele lê todos estes snapshots e permite navegar no tempo.

---

## 8. Roadmap dos outputs

| Output | Fase | Estado |
|---|---|---|
| U (velocidade) | A | Disponível |
| p (pressão) | A | Disponível |
| WSS (Wall Shear Stress) | C | Planeado |
| TAWSS (Time-Averaged WSS) | C | Planeado |
| OSI (Oscillatory Shear Index) | C | Planeado |
| Gradiente de WSS | C | Planeado |
| Relatório JSON com métricas | C | Planeado |
| Screenshots ParaView automáticos | C | Planeado |
| Carreau (viscosidade não-Newtoniana) | D | Planeado |
