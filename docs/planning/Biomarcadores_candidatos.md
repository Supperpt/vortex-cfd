# Biomarcadores Candidatos para a Predição de Rutura em SIAs (<5mm)
## Eixo Hemodinâmico: O Paradigma do "Alto Fluxo" e Instabilidade
Ao contrário dos aneurismas grandes, onde o baixo WSS e a inflamação predominam, os SIAs parecem romper devido ao estresse mecânico direto e instabilidade direcional.

- **Oscillatory Shear Index (OSI):**
	Por que é provável: Foi identificado como um discriminante independente e específico para aneurismas <5mm, enquanto o baixo WSS foi mais relevante para os grandes.
	Evidência: Um aumento unitário no OSI pode aumentar as chances de rutura em até 17.4 vezes. Pequenos aneurismas rotos apresentam áreas de OSI elevado localizadas exatamente no ponto de rutura.
- **High Shear Concentration Ratio (HSCR) e WSS Máximo:**
	Por que é provável: Captura a hipótese de "excesso de energia". O HSCR foca na concentração local da força, sendo significativamente superior no grupo que apresenta instabilidade em aneurismas de 3–5 mm.
	Evidência: Regiões de parede fina (TWRs) — os pontos de rutura iminente — exibem WSS e TaWSS significativamente mais elevados (+3.3%) do que as áreas normais do saco.
- **Wall Shear Stress Divergence (WSSD):**
	Por que é provável: É considerado o indicador mais robusto de instabilidade e adelgaçamento da parede em aneurismas pequenos, capturando as forças de tração do impacto do fluxo melhor que o WSS médio.
- **Oscillatory Velocity Index (OVI):**
	Por que é provável: Quantifica a instabilidade do fluxo em 3D. Em estudos de pares, o OVI mostrou-se significativamente superior em aneurismas rotos mesmo quando as dimensões primárias eram idênticas.

## Eixo Morfológico: Superando o Tamanho Absoluto
As fontes são unânimes: o diâmetro isolado falha em SIAs, mas a relação com o vaso e a irregularidade da superfície são preditores fortes.

- **Size Ratio (SR):**
	Por que é provável: É consistentemente apontado como o preditor morfológico mais estável e superior ao tamanho absoluto para aneurismas pequenos.
	Evidência: Kashiwazaki et al. (2013) estabeleceram que o SR pode prever rutura em SIAs com alta precisão, definindo um limiar de corte de 3.12.
- **Non-Sphericity Index (NSI) e Undulation Index (UI):**
	Por que é provável: Medem o desvio de uma forma esférica e a irregularidade da superfície. O NSI é frequentemente o melhor preditor de forma (AUC ~0.83), indicando que quanto mais o SIA "estica" ou se torna irregular, maior o risco.
	Nota: Em alguns modelos de 2025 para SIAs, o NSI e o Neck Inflow Rate foram os únicos preditores hemodinâmicos/morfológicos a atingir significância inicial.

## Padrões de Fluxo Qualitativos

- **Flow Complexity e Stability:**
	Aneurismas rotos (<5mm) têm uma probabilidade muito maior de apresentar padrões de fluxo complexos (múltiplos vórtices) e fluxo instável durante o ciclo cardíaco.
	Flow Concentration Ratio (FCR): A presença de um jato de entrada concentrado que atinge uma área pequena da parede (zona de impacto) é um marcador clássico de rutura em pequenas dimensões.

## Resumo de Probabilidades para a sua Amostra (<5mm)

| Biomarcador | Tipo | Probabilidade de Significância | Justificação |
|---|---|---|---|
| Size Ratio (SR) | Morfológico | Altíssima | Integra a influência do vaso progenitor no risco local |
| OSI | Hemodinâmico | Altíssima | Principal discriminante específico para o tamanho <5mm |
| HSCR / WSSD | Hemodinâmico | Alta | Identifica a falha da parede por excesso de carga pontual.
| NSI | Morfológico | Alta | Captura a perda de estabilidade geométrica do saco.
| Complex Flow | Hemodinâmico | Média/Alta | Vórtices múltiplos e instáveis são típicos de casos rotos.


Nota Crítica de Cautela (O "Aviso Swiatek"):
Deve ter em atenção que, especificamente em aneurismas <5mm, existe um risco elevado de overfitting. Modelos que mostram excelente desempenho em conjuntos de treino podem colapsar em validações externas (ex: queda de pseudo-R² de 0.142 para 0.032). Para que os seus resultados sejam robustos, a validação em dados independentes é tão importante quanto a escolha do biomarcador.
