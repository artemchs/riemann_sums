# Calculating Riemann Sums

In the following exercises, we graph each function over the given interval, partition the interval into four subintervals of equal length, and add the rectangles associated with the Riemann sum — given that the sample point is the **(a) left-hand endpoint**, **(b) right-hand endpoint**, **(c) midpoint** of the k-th subinterval.

> The full report with plots is available in [`calculating_sums.pdf`](./calculating_sums.pdf).

## Setup

**Constant:** `n = 4` subintervals

### Helper functions

```julia
using CairoMakie

n = 4

function riemann_plot(f, a, b, cpoints; title="")
    n = length(cpoints)
    Δ = (b - a) / n
    xs = range(a, b, length=800)

    fig = Figure(size=(900, 500))
    ax = Axis(fig[1, 1], title=title, xlabel="x", ylabel="y")
    hlines!(ax, [0], color=:gray70, linewidth=1)

    lines!(ax, xs, f.(xs), color=:black, linewidth=3)

    for (k, c) in enumerate(cpoints)
        left = a + (k - 1) * Δ
        right = a + k * Δ
        h = f(c)

        poly!(
            ax,
            Point2f[(left, 0), (left, h), (right, h), (right, 0)],
            color=RGBAf(0.2, 0.5, 0.9, 0.22),
            strokecolor=:dodgerblue,
            strokewidth=2
        )

        scatter!(ax, [c], [h], color=:crimson, markersize=12)
        vlines!(ax, [c], color=(:crimson, 0.25), linestyle=:dash)
    end

    fig
end

partition(a, Δ) = [a + k * Δ for k in 0:n]
left_hand_endpoint(p) = [p[k] for k in 1:n]
right_hand_endpoint(p) = [p[k] for k in 2:(n + 1)]
midpoint(p, Δ) = [p[k] + (Δ / 2) for k in 1:n]
```

---

## Exercise 1 — f(x) = x² − 1 on [0, 2]

```julia
f1(x) = x^2 - 1
a1, b1 = 0, 2
Δ1 = (b1 - a1) / n          # 0.5
partition1 = partition(a1, Δ1)  # [0.0, 0.5, 1.0, 1.5, 2.0]
```

**(a) Left-hand endpoint** — sample points: `[0.0, 0.5, 1.0, 1.5]`
```julia
ex1a_cpoints = left_hand_endpoint(partition1)
riemann_plot(f1, a1, b1, ex1a_cpoints; title="Ex. 1: a")
```

**(b) Right-hand endpoint** — sample points: `[0.5, 1.0, 1.5, 2.0]`
```julia
ex1b_cpoints = right_hand_endpoint(partition1)
riemann_plot(f1, a1, b1, ex1b_cpoints; title="Ex. 1: b")
```

**(c) Midpoint** — sample points: `[0.25, 0.75, 1.25, 1.75]`
```julia
ex1c_cpoints = midpoint(partition1, Δ1)
riemann_plot(f1, a1, b1, ex1c_cpoints; title="Ex. 1: c")
```

---

## Exercise 2 — f(x) = −x² on [0, 1]

```julia
f2(x) = -x^2
a2, b2 = 0, 1
Δ2 = (b2 - a2) / n          # 0.25
partition2 = partition(a2, Δ2)  # [0.0, 0.25, 0.5, 0.75, 1.0]
```

**(a) Left-hand endpoint** — sample points: `[0.0, 0.25, 0.5, 0.75]`
```julia
ex2a_cpoints = left_hand_endpoint(partition2)
riemann_plot(f2, a2, b2, ex2a_cpoints; title="Ex. 2: a")
```

**(b) Right-hand endpoint** — sample points: `[0.25, 0.5, 0.75, 1.0]`
```julia
ex2b_cpoints = right_hand_endpoint(partition2)
riemann_plot(f2, a2, b2, ex2b_cpoints; title="Ex. 2: b")
```

**(c) Midpoint** — sample points: `[0.125, 0.375, 0.625, 0.875]`
```julia
ex2c_cpoints = midpoint(partition2, Δ2)
riemann_plot(f2, a2, b2, ex2c_cpoints; title="Ex. 2: c")
```

---

## Exercise 3 — f(x) = sin(x) on [−π, π]

```julia
f3(x) = sin(x)
a3, b3 = -π, π
Δ3 = (b3 - a3) / n          # ≈ 1.5708
partition3 = partition(a3, Δ3)  # [-π, -π/2, 0, π/2, π]
```

**(a) Left-hand endpoint** — sample points: `[-π, -π/2, 0.0, π/2]`
```julia
ex3a_cpoints = left_hand_endpoint(partition3)
riemann_plot(f3, a3, b3, ex3a_cpoints; title="Ex. 3: a")
```

**(b) Right-hand endpoint** — sample points: `[-π/2, 0.0, π/2, π]`
```julia
ex3b_cpoints = right_hand_endpoint(partition3)
riemann_plot(f3, a3, b3, ex3b_cpoints; title="Ex. 3: b")
```

**(c) Midpoint** — sample points: `[-3π/4, -π/4, π/4, 3π/4]`
```julia
ex3c_cpoints = midpoint(partition3, Δ3)
riemann_plot(f3, a3, b3, ex3c_cpoints; title="Ex. 3: c")
```

---

## Exercise 4 — f(x) = sin(x) + 1 on [−π, π]

```julia
f4(x) = sin(x) + 1
a4, b4 = -π, π
Δ4 = (b4 - a4) / n          # ≈ 1.5708
partition4 = partition(a4, Δ4)  # [-π, -π/2, 0, π/2, π]
```

**(a) Left-hand endpoint** — sample points: `[-π, -π/2, 0.0, π/2]`
```julia
ex4a_cpoints = left_hand_endpoint(partition4)
riemann_plot(f4, a4, b4, ex4a_cpoints; title="Ex. 4: a")
```

**(b) Right-hand endpoint** — sample points: `[-π/2, 0.0, π/2, π]`
```julia
ex4b_cpoints = right_hand_endpoint(partition4)
riemann_plot(f4, a4, b4, ex4b_cpoints; title="Ex. 4: b")
```

**(c) Midpoint** — sample points: `[-3π/4, -π/4, π/4, 3π/4]`
```julia
ex4c_cpoints = midpoint(partition4, Δ4)
riemann_plot(f4, a4, b4, ex4c_cpoints; title="Ex. 4: c")
```

---

## Finding the Norm of a Partition

The **norm** of a partition P, written ‖P‖, is the largest of all the subinterval widths.

### Example 5

```julia
P1 = [0, 1.2, 1.5, 2.3, 2.6, 3]
P1_widths = diff(P1)   # [1.2, 0.3, 0.8, 0.3, 0.4]
P1_norm = maximum(P1_widths)  # 1.2
```

### Example 6

```julia
P2 = [-2, -1.6, -0.5, 0, 0.8, 1]
P2_widths = diff(P2)   # [0.4, 1.1, 0.5, 0.8, 0.2]
P2_norm = maximum(P2_widths)  # 1.1
```
