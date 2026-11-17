# ds000008 sub-01: agreement with `main`

AR(1), all 380 trials: the exact method had maximum double absolute beta
difference 1.42e-12 before the float32 write; no trial exceeded the 1e-10
double tolerance. Other parent modes were compared using saved float32 beta
images, 380 trials per mode:

| Mode / method | Bitwise-equal files | Maximum difference | Voxel entries differing (of ~83 M) | NaN pattern |
|---|---:|---|---:|---|
| None / exact | 191/380 | 4.8e-7, 1 float32 ULP | 252 | same |
| FAST / exact | 114/380 | 1 float32 ULP | 568 | same |
| rWLS / exact | 380/380 | 0 | 0 | same |

These compare saved float32 outputs. A separate end-to-end
AR(1) check gave 22 ULP at maximum against main and passed its relative L2
criterion; see [end-to-end checks](integration.md).
