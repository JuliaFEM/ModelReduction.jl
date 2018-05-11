# Theory

The Craig-Bampton method is a dynamic reduction technique that reduces
the mass and stiffness matrices of the model by expressing the boundary
modes in physical coordinates and the elastic modes in modal coordinates.

The equation of motion is:

```math
\begin{equation}
\boldsymbol{M}\ddot{\boldsymbol{u}}+\boldsymbol{K}\boldsymbol{u}=\boldsymbol{f}
\end{equation}
```

The matrices are partitioned into boundary nodes R and the independent
elastic nodes L:

```math
\begin{equation}
\boldsymbol{u}=\begin{bmatrix}\boldsymbol{u}_{\mathrm{R}}\\
\boldsymbol{u}_{\mathrm{L}}
\end{bmatrix}
\end{equation}
```

Equation (1) becomes:

```math
\begin{equation}
\left[\begin{array}{cc}
\boldsymbol{M}_{\mathrm{R}\mathrm{R}} & \boldsymbol{M}_{\mathrm{R}\mathrm{L}}\\
\boldsymbol{M}_{\mathrm{L}\mathrm{R}} & \boldsymbol{M}_{\mathrm{L}\mathrm{L}}
\end{array}\right]\left[\begin{array}{c}
\ddot{\boldsymbol{u}}_{\mathrm{R}}\\
\ddot{\boldsymbol{u}}_{\mathrm{L}}
\end{array}\right]+\left[\begin{array}{cc}
\boldsymbol{K}_{\mathrm{R}\mathrm{R}} & \boldsymbol{K}_{\mathrm{R}\mathrm{L}}\\
\boldsymbol{K}_{\mathrm{L}\mathrm{R}} & \boldsymbol{K}_{\mathrm{L}\mathrm{L}}
\end{array}\right]\left[\begin{array}{c}
\boldsymbol{u}_{\mathrm{R}}\\
\boldsymbol{u}_{\mathrm{L}}
\end{array}\right]=\left[\begin{array}{c}
\boldsymbol{f}_{\mathrm{R}}\\
\boldsymbol{f}_{\mathrm{L}}
\end{array}\right]
\end{equation}
```

The degrees of freedom are are transformed to hybrid coordinates

```math
\begin{equation}
\begin{bmatrix}\boldsymbol{u}_{\mathrm{R}}\\
\boldsymbol{u}_{\mathrm{L}}
\end{bmatrix}=\begin{bmatrix}\boldsymbol{I} & \boldsymbol{0}\\
\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{X}_{\mathrm{L}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{u}_{\mathrm{R}}\\
\boldsymbol{q}_{\mathrm{m}}
\end{bmatrix}
\end{equation}
```

Equation (1) can be rewritten as

```math
\begin{equation}
\begin{bmatrix}\boldsymbol{M}_{\mathrm{RR}} & \boldsymbol{M}_{\mathrm{\mathrm{R}L}}\\
\boldsymbol{M}_{\mathrm{LR}} & \boldsymbol{M}_{\mathrm{LL}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{I} & \boldsymbol{0}\\
\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{X}_{\mathrm{L}}
\end{bmatrix}\begin{bmatrix}\ddot{\boldsymbol{u}}_{\mathrm{R}}\\
\ddot{\boldsymbol{q}}_{\mathrm{m}}
\end{bmatrix}+\begin{bmatrix}\boldsymbol{K}_{\mathrm{RR}} & \boldsymbol{K}_{\mathrm{RL}}\\
\boldsymbol{K}_{\mathrm{LR}} & \boldsymbol{K}_{\mathrm{LL}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{I} & \boldsymbol{0}\\
\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{X}_{\mathrm{L}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{u}_{\mathrm{R}}\\
\boldsymbol{q}_{\mathrm{m}}
\end{bmatrix}=\begin{bmatrix}\boldsymbol{f}_{\mathrm{R}}\\
\boldsymbol{0}
\end{bmatrix}
\end{equation}
```
Equation (1) reduces to

```math
\begin{equation}
\boldsymbol{K}_{\mathrm{LR}}\boldsymbol{K}_{\mathrm{LR}}\boldsymbol{u}_{\mathrm{R}}+\boldsymbol{K}_{\mathrm{LL}}\boldsymbol{u}_{\mathrm{L}}
\end{equation}
```

The internal degrees of freedom can be expressed as

```math
\begin{equation}
\boldsymbol{u}_{\mathrm{L}}=-\boldsymbol{K}_{\mathrm{LL}}^{-1}\boldsymbol{K}_{\mathrm{LR}}\boldsymbol{u}_{\mathrm{R}}=\boldsymbol{X}_{\mathrm{R}}\boldsymbol{u}_{\mathrm{R}}
\end{equation}
```

where

```math
\begin{equation}
\boldsymbol{X}_{\mathrm{R}}=-\boldsymbol{K}_{\mathrm{LL}}^{-1}\boldsymbol{K}_{\mathrm{LR}}
\end{equation}
```

To determine $\mathit{\boldsymbol{X}}_{\mathrm{L}}$ the retained degrees of freedom are fixed. The equation of motion reduces to

```math
\begin{equation}
\boldsymbol{M}_{\mathrm{LL}}\ddot{\boldsymbol{u}}_{\mathrm{L}}+\boldsymbol{K}_{\mathrm{LL}}\boldsymbol{u}_{\mathrm{L}}=0
\end{equation}
```

By assuming harmonic response and substituting the coordinate transformation (4)

```math
\begin{equation}
(-\omega^{2}\boldsymbol{M}_{\mathrm{LL}}+\boldsymbol{K}_{\mathrm{LL}})\boldsymbol{X}_{\mathrm{L}}\boldsymbol{q}_{\mathrm{m}}e^{i\omega t}=0
\end{equation}
```

The eigenvectors can be normalized:

```math
\begin{equation}
\boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}\boldsymbol{M}_{\mathrm{LL}}\boldsymbol{X}_{\mathrm{L}}=\boldsymbol{I}
\end{equation}
```

```math
\begin{equation}
\boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}\boldsymbol{K}_{\mathrm{LL}}\boldsymbol{X}_{\mathrm{L}}=\boldsymbol{\Lambda}
\end{equation}
```

Since $\boldsymbol{X}$$_{\mathrm{R}}$ in (9) contains $\boldsymbol{K}$$_{\mathrm{LL}}^{-1}$,
an inverse of $\boldsymbol{K}$$_{\mathrm{LL}}$, determining it will
require lots of computing resources. This can be avoided by determining
the $\boldsymbol{K}$$_{\mathrm{LL}}$ inverse as follows.

```math
\begin{equation}
-\boldsymbol{K}_{\mathrm{LL}}^{-1}=\boldsymbol{X}_{\mathrm{L}}\boldsymbol{\Lambda}^{-1}\boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}
\end{equation}
```

In order to get the dynamic equations of the system, equation (6) is multiplied with the coordination transformation matrix.

```math
\begin{multline}
\begin{bmatrix}\boldsymbol{I} & \boldsymbol{X}_{\mathrm{R}}^{\mathrm{T}}\\
\boldsymbol{0} & \boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{M}_{\mathrm{RR}} & \boldsymbol{M}_{\mathrm{\mathrm{R}L}}\\
\boldsymbol{M}_{\mathrm{LR}} & \boldsymbol{M}_{\mathrm{LL}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{I} & \boldsymbol{0}\\
\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{X}_{\mathrm{L}}
\end{bmatrix}\begin{bmatrix}\ddot{\boldsymbol{u}}_{\mathrm{R}}\\
\ddot{\boldsymbol{q}}_{\mathrm{m}}
\end{bmatrix}\\
+\begin{bmatrix}\boldsymbol{I} & \boldsymbol{X}_{\mathrm{R}}^{\mathrm{T}}\\
\boldsymbol{0} & \boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{K}_{\mathrm{RR}} & \boldsymbol{K}_{\mathrm{RL}}\\
\boldsymbol{K}_{\mathrm{LR}} & \boldsymbol{K}_{\mathrm{LL}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{I} & \boldsymbol{0}\\
\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{X}_{\mathrm{L}}
\end{bmatrix}\begin{bmatrix}\mathbf{\mathit{\boldsymbol{u}}}_{\mathrm{R}}\\
\boldsymbol{q}_{\mathrm{m}}
\end{bmatrix}=\begin{bmatrix}\boldsymbol{I} & \boldsymbol{X}_{\mathrm{R}}^{\mathrm{T}}\\
\boldsymbol{0} & \boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}
\end{bmatrix}\begin{bmatrix}\boldsymbol{f}_{\mathrm{R}}\\
\boldsymbol{0}
\end{bmatrix}
\end{multline}
```

By simplifying the equation of motion (1) becomes

```math
\begin{multline}
\begin{bmatrix}\boldsymbol{M}_{\mathrm{RR}}+\mathbf{M}_{\mathbf{RL}}\mathbf{X}_{\mathbf{R}}+\boldsymbol{X}_{\mathrm{R}}^{\mathrm{T}}\boldsymbol{M}_{\mathrm{LR}}+\boldsymbol{X}_{\mathrm{R}}^{\mathrm{T}}\boldsymbol{M}_{\mathrm{LL}}\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{M}_{\mathrm{\mathrm{R}L}}\boldsymbol{X}_{\mathrm{L}}+\boldsymbol{X}_{\mathrm{R}}^{\mathrm{T}}\boldsymbol{M}_{\mathrm{LL}}\boldsymbol{X}_{\mathrm{L}}\\
\boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}\boldsymbol{M}_{\mathrm{LR}}+\boldsymbol{X}_{\mathrm{L}}^{\mathrm{T}}\boldsymbol{M}_{\mathrm{LL}}\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{I}
\end{bmatrix}\begin{bmatrix}\ddot{\boldsymbol{u}}_{\mathrm{R}}\\
\ddot{\boldsymbol{q}}_{\mathrm{m}}
\end{bmatrix}\\
+\begin{bmatrix}\boldsymbol{K}_{\mathrm{RR}}+\boldsymbol{K}_{\mathrm{RL}}\boldsymbol{X}_{\mathrm{R}} & \boldsymbol{0}\\
\boldsymbol{0} & \boldsymbol{\Lambda}
\end{bmatrix}\begin{bmatrix}\boldsymbol{u}_{\mathrm{R}}\\
\boldsymbol{q}_{\mathrm{m}}
\end{bmatrix}=\begin{bmatrix}\boldsymbol{f}_{\mathrm{R}}\\
\boldsymbol{0}
\end{bmatrix}
\end{multline}
```

# References

- Qu, Zu-Qing. Model Order Reduction Techniques (2004). p. 322 - 329.
- Haile, William B. Prime on the Craig-Bampton Method (2000). p. 5 - 17.
