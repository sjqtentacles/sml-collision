# sml-collision

2D AABB, circle, and SAT convex polygon collision detection in pure Standard ML.

Built on [sml-glm](https://github.com/sjqtentacles/sml-glm) for vector math.
Pure and total: no global state, no exceptions for ordinary inputs. Builds with
both [MLton](http://mlton.org/) and [Poly/ML](https://polyml.org/).

## Features

- **AABB** vs AABB — interval overlap on both axes.
- **Circle** vs circle — squared-distance test against the sum of radii.
- **AABB vs circle** — clamp the circle center to the box, compare distance to radius.
- **Convex polygons (CCW)** via the **Separating Axis Theorem (SAT)**, including
  polygon-vs-AABB (the AABB is treated as a 4-point polygon).
- **Segments** as degenerate shapes, with point/segment-vs-box and -vs-circle tests.
- **Collision manifolds** (`collide`): contact normal, positive penetration depth,
  and an approximate contact point.
- **Spatial hash `Grid`** — a pure, persistent broad-phase that buckets shapes by
  cell and returns candidate ids overlapping a query shape's AABB.

Boundary contact is treated as overlapping (closed tests): two AABBs sharing an
edge, or two circles exactly `r1 + r2` apart, both report `true`.

## Installation

```
smlpkg add github.com/sjqtentacles/sml-collision
smlpkg sync
```

The vendored `sml-glm` dependency is included under
`lib/github.com/sjqtentacles/sml-glm`.

## Usage

```sml
open Collision
structure V = Collision.Glm.Vec2

(* AABB overlap *)
val a = AABB { min = V.v (0.0, 0.0), max = V.v (1.0, 1.0) }
val b = AABB { min = V.v (0.5, 0.5), max = V.v (1.5, 1.5) }
val hit = overlap a b                       (* true *)

(* Penetration manifold *)
val SOME { normal, depth, contacts } =
      collide a (AABB { min = V.v (0.5, 0.0), max = V.v (1.5, 1.0) })
(* depth ~ 0.5, normal along +x *)

(* Circles *)
val c1 = Circle { center = V.v (0.0, 0.0), radius = 1.0 }
val c2 = Circle { center = V.v (1.5, 0.0), radius = 1.0 }
val touching = overlap c1 c2                 (* true *)

(* Convex polygons via SAT (vertices CCW) *)
val tri = Poly [ V.v (0.0, 0.0), V.v (2.0, 0.0), V.v (1.0, 1.7) ]
val selfHit = overlap tri tri                (* true *)

(* Broad-phase spatial hash *)
val g = Grid.insert (Grid.insert (Grid.make 1.0) 1 a) 2
          (AABB { min = V.v (10.0, 10.0), max = V.v (11.0, 11.0) })
val candidates = Grid.query g a              (* [1] *)
```

## API

```sml
type vec2 = Glm.Vec2.t
type aabb    = { min : vec2, max : vec2 }
type circle  = { center : vec2, radius : real }
type segment = { a : vec2, b : vec2 }
datatype shape = AABB of aabb | Circle of circle | Segment of segment | Poly of vec2 list
type manifold = { normal : vec2, depth : real, contacts : vec2 list }

val overlap : shape -> shape -> bool
val collide : shape -> shape -> manifold option

structure Grid :
sig
  type t
  val make   : real -> t          (* cell size *)
  val insert : t -> int -> shape -> t
  val remove : t -> int -> t
  val query  : t -> shape -> int list
  val clear  : t -> t
end
```

## Testing

```
make test       # MLton
make test-poly  # Poly/ML
make all-tests  # both
```

## License

MIT
