    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Masters that need to be restored:56
Edge intersection testing:
    Number of edges             : 1624999
    Number of edges to retest   : 394
    Number of intersected edges : 134784


Undo iteration 1
----------------
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0

Merging all points on surface that
- are used by only two boundary faces and
- make an angle with a cosine of more than 0.5.

Removing 1220 straight edge points ...

Edge intersection testing:
    Number of edges             : 1624999
    Number of edges to retest   : 8977
    Number of intersected edges : 134780

Undo iteration 0
----------------
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 1
    faces with face-decomposition tet quality < 1e-30      : 4
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Detected 4 error faces on boundaries that have been merged. These will be restored to their original faces.

Edge intersection testing:
    Number of edges             : 1624999
    Number of edges to retest   : 34
    Number of intersected edges : 134780

Undo iteration 1
----------------
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 1
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Detected 0 error faces on boundaries that have been merged. These will be restored to their original faces.

Detected 1 error faces in mesh. Restoring neighbours of faces in error.

Edge intersection testing:
    Number of edges             : 1624999
    Number of edges to retest   : 62
    Number of intersected edges : 134780

Checking mesh manifoldness ...

Checking initial mesh ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Detected 0 illegal faces (concave, zero area or negative cell pyramid volume)
Duplicating 0 points on faceZones of type boundary

Handling points with inconsistent layer specification ...
dupNonManifoldPoints : Found : 0 non-manifold points (out of 651027)

Handling points with inconsistent layer specification ...

Handling non-manifold points ...

Checking patch manifoldness ...
Outside of local patch is multiply connected across edges or points at 0 points.
Set displacement to zero for all 0 non-manifold points

Handling feature edges (angle < 60) ...
Set displacement to zero for points on 0 feature edges

Handling cells with warped patch faces ...
Set displacement to zero on 0 warped faces since layer would be > 0.5 of the size of the bounding box.

patch         faces    layers avg thickness[m]
                              near-wall overall
-----         -----    ------ --------- -------
aneurysm_sac  8432     4      6e-05     0.000371
parent_vessel 122458   4      6.65e-05  0.000412


Outer iteration : 0
-------------------

Adding in total 0 inter-processor patches to handle extrusion of non-manifold processor boundaries.
Selecting externalDisplacementMeshMover displacementMedialAxis
displacementMedialAxis : Calculating distance to Medial Axis ...
fieldSmoother : Smoothing normals ...
    Iteration 0   residual 0.010430256
displacementMedialAxis : Inserting points on patch background if angle to nearest layer patch > 30 degrees.
displacementMedialAxis : Inserting points on patch inlet if angle to nearest layer patch > 30 degrees.
displacementMedialAxis : Inserting points on patch outlet_0 if angle to nearest layer patch > 30 degrees.
displacementMedialAxis : Inserting points on patch outlet_1 if angle to nearest layer patch > 30 degrees.
fieldSmoother : Smoothing normals in interior ...
    Iteration 0   residual 0.023082715

Layer addition iteration 0
--------------------------

Determining displacement for added points according to pointNormal ...
Detected 0 points with point normal pointing through faces.
Reset displacement at 0 points to average of surrounding points.

displacementMedialAxis : Smoothing using Medial Axis ...
displacementMedialAxis : Reducing layer thickness at 0 nodes where thickness to medial axis distance is large 
displacementMedialAxis : Removing isolated regions ...
- if partially extruded faces make angle < 30
- if exclusively surrounded by non-extruded faces
displacementMedialAxis : Number of isolated points extrusion stopped : 0
fieldSmoother : Smoothing field ...
    Iteration 0   residual 4.1435219e-07
displacementMedialAxis : Moving mesh ...
displacementMedialAxis : Iteration 0
Moving mesh using displacement scaling : min:1  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 153
    faces with concavity > 80  degrees                     : 2
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 1
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 1
Moving mesh using displacement scaling : min:0.75  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 99
    faces with concavity > 80  degrees                     : 1
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 1
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 2
Moving mesh using displacement scaling : min:0.5625  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 65
    faces with concavity > 80  degrees                     : 1
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 1
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 3
displacementMedialAxis : Displacement scaling for error reduction set to 0.
Moving mesh using displacement scaling : min:0.421875  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 56
    faces with concavity > 80  degrees                     : 1
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 4
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 13
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 5
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 6
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Finished moving mesh ...
truncateDisplacement : Unextruded 1 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.
truncateDisplacement : Unextruded 0 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.

Setting up information for layer truncation ...
Detected 0 baffles across faceZones of type internal


Checking mesh with layer ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 26
    faces with face pyramid volume < 1e-13                 : 288
    faces with face-decomposition tet quality < 1e-30      : 2255
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 1
    faces with face twist < 0.02                           : 1
    faces on cells with determinant < 0.001                : 1
Detected 2572 illegal faces (concave, zero area or negative cell pyramid volume)
Extruding 129187 out of 130890 faces (98.698907%). Removed extrusion at 1439 faces.
Added 518901 out of 523560 cells (99.110131%).


Layer addition iteration 1
--------------------------

Determining displacement for added points according to pointNormal ...
Detected 0 points with point normal pointing through faces.
Reset displacement at 0 points to average of surrounding points.

displacementMedialAxis : Smoothing using Medial Axis ...
displacementMedialAxis : Reducing layer thickness at 0 nodes where thickness to medial axis distance is large 
displacementMedialAxis : Removing isolated regions ...
- if partially extruded faces make angle < 30
- if exclusively surrounded by non-extruded faces
displacementMedialAxis : Number of isolated points extrusion stopped : 775
fieldSmoother : Smoothing field ...
    Iteration 0   residual 2.370701e-06
displacementMedialAxis : Moving mesh ...
displacementMedialAxis : Iteration 0
Moving mesh using displacement scaling : min:1  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 1
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 90
    faces with concavity > 80  degrees                     : 1
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 1
Moving mesh using displacement scaling : min:0.75  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 1
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 38
    faces with concavity > 80  degrees                     : 1
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 2
Moving mesh using displacement scaling : min:0.5625  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 1
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 29
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 3
displacementMedialAxis : Displacement scaling for error reduction set to 0.
Moving mesh using displacement scaling : min:0.421875  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 1
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 18
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 4
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 2
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 5
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Successfully moved mesh
displacementMedialAxis : Finished moving mesh ...
truncateDisplacement : Unextruded 50 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 3 faces due to stringed edges with inconsistent extrusion.
truncateDisplacement : Unextruded 0 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.

Setting up information for layer truncation ...
Detected 0 baffles across faceZones of type internal


Checking mesh with layer ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 6
    faces with face pyramid volume < 1e-13                 : 62
    faces with face-decomposition tet quality < 1e-30      : 382
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 8
    faces on cells with determinant < 0.001                : 4
Detected 462 illegal faces (concave, zero area or negative cell pyramid volume)
Extruding 128079 out of 130890 faces (97.852395%). Removed extrusion at 408 faces.
Added 462186 out of 523560 cells (88.277561%).


Layer addition iteration 2
--------------------------

Determining displacement for added points according to pointNormal ...
Detected 0 points with point normal pointing through faces.
Reset displacement at 0 points to average of surrounding points.

displacementMedialAxis : Smoothing using Medial Axis ...
displacementMedialAxis : Reducing layer thickness at 0 nodes where thickness to medial axis distance is large 
displacementMedialAxis : Removing isolated regions ...
- if partially extruded faces make angle < 30
- if exclusively surrounded by non-extruded faces
displacementMedialAxis : Number of isolated points extrusion stopped : 376
fieldSmoother : Smoothing field ...
    Iteration 0   residual 2.4285658e-06
displacementMedialAxis : Moving mesh ...
displacementMedialAxis : Iteration 0
Moving mesh using displacement scaling : min:1  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 1
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 67
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 1
Moving mesh using displacement scaling : min:0.75  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 1
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 23
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 2
Moving mesh using displacement scaling : min:0.5625  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 1
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 15
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 3
displacementMedialAxis : Displacement scaling for error reduction set to 0.
Moving mesh using displacement scaling : min:0.421875  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 9
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 4
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Successfully moved mesh
displacementMedialAxis : Finished moving mesh ...
truncateDisplacement : Unextruded 10 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.
truncateDisplacement : Unextruded 0 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.

Setting up information for layer truncation ...
Detected 0 baffles across faceZones of type internal


Checking mesh with layer ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 5
    faces with face pyramid volume < 1e-13                 : 15
    faces with face-decomposition tet quality < 1e-30      : 156
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 3
    faces with face twist < 0.02                           : 5
    faces on cells with determinant < 0.001                : 5
Detected 189 illegal faces (concave, zero area or negative cell pyramid volume)
Extruding 127486 out of 130890 faces (97.399343%). Removed extrusion at 178 faces.
Added 458325 out of 523560 cells (87.54011%).


Layer addition iteration 3
--------------------------

Determining displacement for added points according to pointNormal ...
Detected 0 points with point normal pointing through faces.
Reset displacement at 0 points to average of surrounding points.

displacementMedialAxis : Smoothing using Medial Axis ...
displacementMedialAxis : Reducing layer thickness at 0 nodes where thickness to medial axis distance is large 
displacementMedialAxis : Removing isolated regions ...
- if partially extruded faces make angle < 30
- if exclusively surrounded by non-extruded faces
displacementMedialAxis : Number of isolated points extrusion stopped : 140
fieldSmoother : Smoothing field ...
    Iteration 0   residual 2.413692e-06
displacementMedialAxis : Moving mesh ...
displacementMedialAxis : Iteration 0
Moving mesh using displacement scaling : min:1  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 62
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 1
Moving mesh using displacement scaling : min:0.75  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 21
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 2
Moving mesh using displacement scaling : min:0.5625  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 13
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 3
displacementMedialAxis : Displacement scaling for error reduction set to 0.
Moving mesh using displacement scaling : min:0.421875  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 7
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 4
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Successfully moved mesh
displacementMedialAxis : Finished moving mesh ...
truncateDisplacement : Unextruded 6 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.
truncateDisplacement : Unextruded 0 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.

Setting up information for layer truncation ...
Detected 0 baffles across faceZones of type internal


Checking mesh with layer ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 3
    faces with face pyramid volume < 1e-13                 : 3
    faces with face-decomposition tet quality < 1e-30      : 70
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 1
    faces on cells with determinant < 0.001                : 0
Detected 77 illegal faces (concave, zero area or negative cell pyramid volume)
Extruding 127248 out of 130890 faces (97.217511%). Removed extrusion at 74 faces.
Added 457287 out of 523560 cells (87.341852%).


Layer addition iteration 4
--------------------------

Determining displacement for added points according to pointNormal ...
Detected 0 points with point normal pointing through faces.
Reset displacement at 0 points to average of surrounding points.

displacementMedialAxis : Smoothing using Medial Axis ...
displacementMedialAxis : Reducing layer thickness at 0 nodes where thickness to medial axis distance is large 
displacementMedialAxis : Removing isolated regions ...
- if partially extruded faces make angle < 30
- if exclusively surrounded by non-extruded faces
displacementMedialAxis : Number of isolated points extrusion stopped : 47
fieldSmoother : Smoothing field ...
    Iteration 0   residual 2.3989123e-06
displacementMedialAxis : Moving mesh ...
displacementMedialAxis : Iteration 0
Moving mesh using displacement scaling : min:1  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 62
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 1
Moving mesh using displacement scaling : min:0.75  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 23
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 2
Moving mesh using displacement scaling : min:0.5625  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 13
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 3
displacementMedialAxis : Displacement scaling for error reduction set to 0.
Moving mesh using displacement scaling : min:0.421875  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 7
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 4
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Successfully moved mesh
displacementMedialAxis : Finished moving mesh ...
truncateDisplacement : Unextruded 0 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.

Setting up information for layer truncation ...
Detected 0 baffles across faceZones of type internal


Checking mesh with layer ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 5
    faces with face-decomposition tet quality < 1e-30      : 24
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Detected 29 illegal faces (concave, zero area or negative cell pyramid volume)
Extruding 127155 out of 130890 faces (97.146459%). Removed extrusion at 25 faces.
Added 456971 out of 523560 cells (87.281496%).


Layer addition iteration 5
--------------------------

Determining displacement for added points according to pointNormal ...
Detected 0 points with point normal pointing through faces.
Reset displacement at 0 points to average of surrounding points.

displacementMedialAxis : Smoothing using Medial Axis ...
displacementMedialAxis : Reducing layer thickness at 0 nodes where thickness to medial axis distance is large 
displacementMedialAxis : Removing isolated regions ...
- if partially extruded faces make angle < 30
- if exclusively surrounded by non-extruded faces
displacementMedialAxis : Number of isolated points extrusion stopped : 20
fieldSmoother : Smoothing field ...
    Iteration 0   residual 2.3972305e-06
displacementMedialAxis : Moving mesh ...
displacementMedialAxis : Iteration 0
Moving mesh using displacement scaling : min:1  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 62
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 1
Moving mesh using displacement scaling : min:0.75  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 23
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 2
Moving mesh using displacement scaling : min:0.5625  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 13
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 3
displacementMedialAxis : Displacement scaling for error reduction set to 0.
Moving mesh using displacement scaling : min:0.421875  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 7
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 4
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Successfully moved mesh
displacementMedialAxis : Finished moving mesh ...
truncateDisplacement : Unextruded 2 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.
truncateDisplacement : Unextruded 0 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.

Setting up information for layer truncation ...
Detected 0 baffles across faceZones of type internal


Checking mesh with layer ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 2
    faces with face-decomposition tet quality < 1e-30      : 9
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Detected 11 illegal faces (concave, zero area or negative cell pyramid volume)
Extruding 127116 out of 130890 faces (97.116663%). Removed extrusion at 11 faces.
Added 456804 out of 523560 cells (87.249599%).


Layer addition iteration 6
--------------------------

Determining displacement for added points according to pointNormal ...
Detected 0 points with point normal pointing through faces.
Reset displacement at 0 points to average of surrounding points.

displacementMedialAxis : Smoothing using Medial Axis ...
displacementMedialAxis : Reducing layer thickness at 0 nodes where thickness to medial axis distance is large 
displacementMedialAxis : Removing isolated regions ...
- if partially extruded faces make angle < 30
- if exclusively surrounded by non-extruded faces
displacementMedialAxis : Number of isolated points extrusion stopped : 8
fieldSmoother : Smoothing field ...
    Iteration 0   residual 2.3942051e-06
displacementMedialAxis : Moving mesh ...
displacementMedialAxis : Iteration 0
Moving mesh using displacement scaling : min:1  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 62
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 1
Moving mesh using displacement scaling : min:0.75  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 23
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 2
Moving mesh using displacement scaling : min:0.5625  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 13
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 3
displacementMedialAxis : Displacement scaling for error reduction set to 0.
Moving mesh using displacement scaling : min:0.421875  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 7
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Iteration 4
Moving mesh using displacement scaling : min:0  max:1
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
displacementMedialAxis : Successfully moved mesh
displacementMedialAxis : Finished moving mesh ...
truncateDisplacement : Unextruded 0 faces due to non-consecutive vertices being extruded.
truncateDisplacement : Unextruded 0 faces due to stringed edges with inconsistent extrusion.

Setting up information for layer truncation ...
Detected 0 baffles across faceZones of type internal


Checking mesh with layer ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Detected 0 illegal faces (concave, zero area or negative cell pyramid volume)
Extruding 127106 out of 130890 faces (97.109023%). Removed extrusion at 0 faces.
Added 456762 out of 523560 cells (87.241577%).
Edge intersection testing:
    Number of edges             : 3120099
    Number of edges to retest   : 0
    Number of intersected edges : 594198
Mesh with layers : cells:961450  faces:3120099  points:1232091

patch         faces        layers        overall thickness
                       target   mesh     [m]       [%]
-----         -----    -----    ----     ---       ---
aneurysm_sac  8432     4        2.72     0.000293  74.4    
parent_vessel 122458   4        3.54     0.000384  91      

Writing 456762 added cells to cellSet addedCells
Writing 329656 faces inside added layer to faceSet layerFaces

Writing fields with layer information:
    nSurfaceLayers    : actual number of layers
    thickness         : overall layer thickness
    thicknessFraction : overall layer thickness (fraction of desired thickness)

Layer mesh : cells:961450  faces:3120099  points:1232091
Cells per refinement level:
    0   103
    1   9418
    2   67781
    3   812852
    4   71296
Writing mesh to time constant
Wrote mesh in = 0.05 s.
Layers added in = 65.9 s.
Checking final mesh ...
Checking faces in error :
    non-orthogonality > 70  degrees                        : 0
    faces with face pyramid volume < 1e-13                 : 0
    faces with face-decomposition tet quality < 1e-30      : 0
    faces with concavity > 80  degrees                     : 0
    faces with skewness > 4   (internal) or 20  (boundary) : 0
    faces with interpolation weights (0..1)  < 0.02        : 0
    faces with volume ratio of neighbour cells < 0.01      : 0
    faces with face twist < 0.02                           : 0
    faces on cells with determinant < 0.001                : 0
Finished meshing without any errors
Finished meshing in = 206.05 s.
End


[vortex-cfd] checkMesh
  $ checkMesh
/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2406                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
Build  : _1653fa08-20260127 OPENFOAM=2406 patch=260127 version=2406
Arch   : "LSB;label=32;scalar=64"
Exec   : checkMesh
Date   : Jun 13 2026
Time   : 16:05:50
Host   : LUZPC
PID    : 21512
I/O    : uncollated
Case   : /home/diogovluz/Documents/simulations/AA_011/case_20260613_160223
nProcs : 1
trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 5, maxFileModificationPolls 20)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

Create mesh for time = 0

Check mesh...

Time = 0

Mesh stats 
    points:           1232091
    faces:            3120099
    internal faces:   2982562
    cells:            961450
    faces per cell:   6.3473514
    boundary patches: 5
    point zones:      1
    face zones:       0
    cell zones:       0

Overall number of cells of each type:
    hexahedra:     662153
    prisms:        17303
    wedges:        4
    pyramids:      0
    tet wedges:    707
    tetrahedra:    9
    polyhedra:     281274
    Breakdown of polyhedra by number of faces:
        faces   number of cells
            4   34752
            5   15684
            6   19704
            7   112702
            8   57285
            9   19863
           10   2456
           11   670
           12   8960
           13   44
           14   50
           15   8626
           16   9
           17   1
           18   458
           21   10

Checking topology...
    Boundary definition OK.
    Cell to face addressing OK.
    Point usage OK.
    Upper triangular ordering OK.
    Face vertices OK.
    Number of regions: 1 (OK).

Checking patch topology for multiply connected surfaces...
                   Patch    Faces   Points  Surface topology
            aneurysm_sac     8432    10693  ok (non-closed singly connected)
                   inlet     2369     3093  ok (non-closed singly connected)
                outlet_0     1074     1446  ok (non-closed singly connected)
                outlet_1     3204     3693  ok (non-closed singly connected)
           parent_vessel   122458   154876  ok (non-closed singly connected)
                    ".*"   137537   172598      ok (closed singly connected)


Checking faceZone topology for multiply connected surfaces...
    No faceZones found.

Checking basic cellZone addressing...
    No cellZones found.

Checking basic pointZone addressing...
               PointZone  PointsBoundingBox
            frozenPoints       0(1e+150 1e+150 1e+150) (-1e+150 -1e+150 -1e+150)

Checking geometry...
    Overall domain bounding box (-0.3058192 0.41468601 -0.13763906) (-0.0192998 0.5924413 0.05113796)
    Mesh has 3 geometric (non-empty/wedge) directions (1 1 1)
    Mesh has 3 solution (non-empty) directions (1 1 1)
    Boundary openness (-3.8291143e-16 2.5804644e-16 6.8734543e-17) OK.
    Max cell openness = 5.1531483e-16 OK.
    Max aspect ratio = 29.477274 OK.
    Minimum face area = 1.3869999e-09. Maximum face area = 2.9806854e-05.  Face area magnitudes OK.
    Min volume = 4.4118425e-13. Max volume = 1.2432794e-07.  Total volume = 0.00037902094.  Cell volumes OK.
    Mesh non-orthogonality Max: 69.959493 average: 12.280944
    Non-orthogonality check OK.
    Face pyramids OK.
 ***Max skewness = 6.6618094, 29 highly skew faces detected which may impair the quality of the results
  <<Writing 29 skew faces to set skewFaces
    Coupled point location match (average 0) OK.

Failed 1 mesh checks.

End

ERROR: maxSkewness = 6.66 > 4 threshold.

Mesh quality check FAILED. Inspect the mesh in ParaView before proceeding.
