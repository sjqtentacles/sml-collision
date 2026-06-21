(* collision.sig

   2D collision detection in pure Standard ML, built on the sml-glm Vec2 type.

   Shapes supported: axis-aligned bounding boxes, circles, segments, and convex
   polygons (vertices in counter-clockwise order). Polygon/polygon and
   polygon/AABB tests use the Separating Axis Theorem (SAT).

   Conventions:
   - `overlap a b` is a pure boolean intersection test. Boundary contact (e.g.
     two AABBs that share an edge, two circles exactly r1+r2 apart) counts as
     overlapping: the tests are closed (>= / <=).
   - `collide a b` returns SOME manifold when the shapes intersect, with a
     contact normal, a positive penetration `depth`, and one or more approximate
     contact points. The normal points from the first shape toward the second.
   - All operations are pure and total; no exceptions are raised for ordinary
     inputs. *)

signature COLLISION =
sig
  (* The vendored sml-glm structure is re-exported so consumers can build
     vectors with `Collision.Glm.Vec2.v` without a separate import. *)
  structure Glm : GLM

  type vec2 = Glm.Vec2.t
  type aabb    = { min : vec2, max : vec2 }
  type circle  = { center : vec2, radius : real }
  type segment = { a : vec2, b : vec2 }

  datatype shape
    = AABB    of aabb
    | Circle  of circle
    | Segment of segment
    | Poly    of vec2 list   (* convex, CCW *)

  type manifold = { normal : vec2, depth : real, contacts : vec2 list }

  val overlap : shape -> shape -> bool
  val collide : shape -> shape -> manifold option

  structure Grid :
  sig
    type t
    val make   : real -> t
    val insert : t -> int -> shape -> t
    val remove : t -> int -> t
    val query  : t -> shape -> int list
    val clear  : t -> t
  end
end
