(* collision.sml — implementation of COLLISION on top of sml-glm Vec2. *)

structure Collision :> COLLISION =
struct
  structure Glm = Glm
  structure V = Glm.Vec2

  type vec2 = V.t
  type aabb    = { min : vec2, max : vec2 }
  type circle  = { center : vec2, radius : real }
  type segment = { a : vec2, b : vec2 }

  datatype shape
    = AABB    of aabb
    | Circle  of circle
    | Segment of segment
    | Poly    of vec2 list

  type manifold = { normal : vec2, depth : real, contacts : vec2 list }

  (* ------------------------------------------------------------------ *)
  (* small numeric helpers *)
  fun fmin (a : real, b) = if a < b then a else b
  fun fmax (a : real, b) = if a > b then a else b
  fun clampf (x, lo, hi) = if x < lo then lo else if x > hi then hi else x

  (* ------------------------------------------------------------------ *)
  (* AABB helpers *)

  (* Normalize a box so min <= max on each axis. *)
  fun normBox { min, max } : aabb =
    let val (x0, x1) = (V.x min, V.x max)
        val (y0, y1) = (V.y min, V.y max)
    in { min = V.v (fmin (x0, x1), fmin (y0, y1)),
         max = V.v (fmax (x0, x1), fmax (y0, y1)) }
    end

  (* Axis-aligned bounding box that encloses a shape. *)
  fun aabbOf (AABB b) = normBox b
    | aabbOf (Circle { center, radius }) =
        let val r = Real.abs radius
        in { min = V.v (V.x center - r, V.y center - r),
             max = V.v (V.x center + r, V.y center + r) }
        end
    | aabbOf (Segment { a, b }) =
        { min = V.v (fmin (V.x a, V.x b), fmin (V.y a, V.y b)),
          max = V.v (fmax (V.x a, V.x b), fmax (V.y a, V.y b)) }
    | aabbOf (Poly pts) =
        (case pts of
           [] => { min = V.zero, max = V.zero }
         | p0 :: rest =>
             let fun step (p, (mnx, mny, mxx, mxy)) =
                       (fmin (mnx, V.x p), fmin (mny, V.y p),
                        fmax (mxx, V.x p), fmax (mxy, V.y p))
                 val (mnx, mny, mxx, mxy) =
                       List.foldl step (V.x p0, V.y p0, V.x p0, V.y p0) rest
             in { min = V.v (mnx, mny), max = V.v (mxx, mxy) } end)

  fun aabbOverlap (a : aabb) (b : aabb) =
    let val a = normBox a and b = normBox b
    in V.x (#min a) <= V.x (#max b) andalso V.x (#max a) >= V.x (#min b)
       andalso V.y (#min a) <= V.y (#max b) andalso V.y (#max a) >= V.y (#min b)
    end

  (* ------------------------------------------------------------------ *)
  (* circle helpers *)

  fun closestOnBox ({ min, max } : aabb) p =
    V.v (clampf (V.x p, V.x min, V.x max), clampf (V.y p, V.y min, V.y max))

  fun aabbCircleOverlap (b : aabb) ({ center, radius } : circle) =
    let val b = normBox b
        val cp = closestOnBox b center
        val d2 = V.lengthSq (V.sub (center, cp))
        val r = Real.abs radius
    in d2 <= r * r
    end

  fun circleOverlap (c1 : circle) (c2 : circle) =
    let val d2 = V.lengthSq (V.sub (#center c1, #center c2))
        val rs = Real.abs (#radius c1) + Real.abs (#radius c2)
    in d2 <= rs * rs
    end

  (* ------------------------------------------------------------------ *)
  (* polygon / SAT helpers *)

  (* The polygon representation used internally for SAT: a vertex list. *)
  fun boxToPoly ({ min, max } : aabb) =
    [ min, V.v (V.x max, V.y min), max, V.v (V.x min, V.y max) ]

  fun shapeToPoly (AABB b)  = SOME (boxToPoly (normBox b))
    | shapeToPoly (Poly ps) = SOME ps
    | shapeToPoly (Segment { a, b }) = SOME [a, b]
    | shapeToPoly (Circle _) = NONE

  (* Edge normals (not normalized) for a polygon's edges. *)
  fun edgeNormals pts =
    let val n = List.length pts
        val arr = Vector.fromList pts
        fun edgeNormal i =
          let val p = Vector.sub (arr, i)
              val q = Vector.sub (arr, (i + 1) mod n)
              val ex = V.x q - V.x p
              val ey = V.y q - V.y p
          in V.v (~ey, ex)   (* perpendicular *)
          end
    in if n < 2 then []
       else List.tabulate (n, edgeNormal)
    end

  (* Project polygon vertices onto axis; return (min, max). *)
  fun project (axis, pts) =
    case pts of
      [] => (0.0, 0.0)
    | p0 :: rest =>
        let val d0 = V.dot (axis, p0)
            fun step (p, (lo, hi)) =
              let val d = V.dot (axis, p) in (fmin (lo, d), fmax (hi, d)) end
        in List.foldl step (d0, d0) rest end

  (* SAT for two convex polygons. Returns SOME (normal, depth) when they
     overlap, where normal is a unit separating axis pointing from a to b
     and depth is the minimum penetration. NONE when a gap exists. *)
  fun satPolys (ptsA, ptsB) =
    let
      val axes = edgeNormals ptsA @ edgeNormals ptsB
      (* centroid direction, used to orient the normal from a toward b *)
      fun centroid pts =
        let val n = Real.fromInt (List.length pts)
            val sx = List.foldl (fn (p, s) => s + V.x p) 0.0 pts
            val sy = List.foldl (fn (p, s) => s + V.y p) 0.0 pts
        in V.v (sx / n, sy / n) end
      val dir = V.sub (centroid ptsB, centroid ptsA)
    in
      (* Any separating axis aborts immediately to NONE. *)
      let
        fun loop2 ([], best) = best
          | loop2 (axis :: rest, best) =
              let val len = V.length axis
              in if len < 1.0E~12 then loop2 (rest, best)
                 else
                   let val n = V.scale (1.0 / len, axis)
                       val (minA, maxA) = project (n, ptsA)
                       val (minB, maxB) = project (n, ptsB)
                   in if maxA < minB orelse maxB < minA
                      then NONE
                      else
                        let val overlapDepth = fmin (maxA, maxB) - fmax (minA, minB)
                            val n = if V.dot (n, dir) < 0.0 then V.neg n else n
                            val better =
                              case best of
                                NONE => SOME (n, overlapDepth)
                              | SOME (_, d) =>
                                  if overlapDepth < d then SOME (n, overlapDepth) else best
                        in loop2 (rest, better) end
                   end
              end
      in loop2 (axes, NONE) end
    end

  (* ------------------------------------------------------------------ *)
  (* segment / point helpers *)

  fun pointInBox (b : aabb) p =
    let val b = normBox b
    in V.x p >= V.x (#min b) andalso V.x p <= V.x (#max b)
       andalso V.y p >= V.y (#min b) andalso V.y p <= V.y (#max b)
    end

  fun pointInCircle ({ center, radius } : circle) p =
    V.lengthSq (V.sub (p, center)) <= Real.abs radius * Real.abs radius

  (* Closest point on segment [a,b] to point p. *)
  fun closestOnSeg (a, b) p =
    let val ab = V.sub (b, a)
        val len2 = V.lengthSq ab
    in if len2 < 1.0E~20 then a
       else
         let val t = clampf (V.dot (V.sub (p, a), ab) / len2, 0.0, 1.0)
         in V.add (a, V.scale (t, ab)) end
    end

  fun segCircleOverlap (a, b) (c : circle) =
    let val cp = closestOnSeg (a, b) (#center c)
    in V.lengthSq (V.sub (#center c, cp)) <= Real.abs (#radius c) * Real.abs (#radius c)
    end

  (* Segment vs AABB: true if either endpoint is inside, or the segment
     crosses any box edge. We treat the box as 4 edges (segments). *)
  fun segSegIntersect (p1, p2) (q1, q2) =
    let fun cross (ax, ay, bx, by) = ax * by - ay * bx
        val r = V.sub (p2, p1)
        val s = V.sub (q2, q1)
        val rxs = cross (V.x r, V.y r, V.x s, V.y s)
        val qp = V.sub (q1, p1)
        val qpxr = cross (V.x qp, V.y qp, V.x r, V.y r)
    in
      if Real.abs rxs < 1.0E~12 then
        (* parallel; treat as non-crossing for our purposes *)
        false
      else
        let val t = cross (V.x qp, V.y qp, V.x s, V.y s) / rxs
            val u = qpxr / rxs
        in t >= 0.0 andalso t <= 1.0 andalso u >= 0.0 andalso u <= 1.0 end
    end

  fun segBoxOverlap (a, b) (bx : aabb) =
    let val bx = normBox bx
    in pointInBox bx a orelse pointInBox bx b
       orelse
       let val c00 = #min bx
           val c10 = V.v (V.x (#max bx), V.y (#min bx))
           val c11 = #max bx
           val c01 = V.v (V.x (#min bx), V.y (#max bx))
       in segSegIntersect (a, b) (c00, c10)
          orelse segSegIntersect (a, b) (c10, c11)
          orelse segSegIntersect (a, b) (c11, c01)
          orelse segSegIntersect (a, b) (c01, c00)
       end
    end

  (* ------------------------------------------------------------------ *)
  (* overlap dispatch *)

  fun overlap s1 s2 =
    case (s1, s2) of
      (AABB a, AABB b) => aabbOverlap (normBox a) (normBox b)
    | (Circle a, Circle b) => circleOverlap a b
    | (AABB b, Circle c) => aabbCircleOverlap b c
    | (Circle c, AABB b) => aabbCircleOverlap b c
    | (Segment { a, b }, Circle c) => segCircleOverlap (a, b) c
    | (Circle c, Segment { a, b }) => segCircleOverlap (a, b) c
    | (Segment { a, b }, AABB bx) => segBoxOverlap (a, b) bx
    | (AABB bx, Segment { a, b }) => segBoxOverlap (a, b) bx
    | (Segment { a = a1, b = b1 }, Segment { a = a2, b = b2 }) =>
        (* degenerate-friendly: coincident endpoints or true crossing *)
        V.equal (a1, a2) orelse V.equal (a1, b2) orelse V.equal (b1, a2)
        orelse V.equal (b1, b2) orelse segSegIntersect (a1, b1) (a2, b2)
    | _ =>
        (* anything involving a Poly (and Poly/Circle) goes through SAT,
           with Circle approximated by sampling its support is overkill;
           for Poly/Circle we fall back to nearest-edge distance. *)
        (case (shapeToPoly s1, shapeToPoly s2) of
           (SOME pa, SOME pb) => Option.isSome (satPolys (pa, pb))
         | (SOME pa, NONE) => polyCircleOverlap pa (asCircle s2)
         | (NONE, SOME pb) => polyCircleOverlap pb (asCircle s1)
         | (NONE, NONE) => circleOverlap (asCircle s1) (asCircle s2))

  and asCircle (Circle c) = c
    | asCircle _ = { center = V.zero, radius = 0.0 }

  (* Polygon vs circle: circle overlaps polygon iff its center is inside the
     polygon or within radius of some edge. *)
  and polyCircleOverlap pts (c : circle) =
    let
      fun pointInPoly p =
        (* winding / half-plane test for convex CCW polygon: inside if on the
           left (>=0) of every edge. We allow either orientation by checking
           consistency of signs. *)
        let val n = List.length pts
            val arr = Vector.fromList pts
            fun side i =
              let val a = Vector.sub (arr, i)
                  val b = Vector.sub (arr, (i + 1) mod n)
                  val e = V.sub (b, a)
                  val w = V.sub (p, a)
              in V.x e * V.y w - V.y e * V.x w end
            val sides = List.tabulate (n, side)
            val allNonNeg = List.all (fn s => s >= ~1.0E~9) sides
            val allNonPos = List.all (fn s => s <= 1.0E~9) sides
        in allNonNeg orelse allNonPos end
      val r2 = Real.abs (#radius c) * Real.abs (#radius c)
      val n = List.length pts
      val arr = Vector.fromList pts
      fun edgeNear i =
        let val a = Vector.sub (arr, i)
            val b = Vector.sub (arr, (i + 1) mod n)
            val cp = closestOnSeg (a, b) (#center c)
        in V.lengthSq (V.sub (#center c, cp)) <= r2 end
    in
      if n = 0 then false
      else pointInPoly (#center c)
           orelse List.exists edgeNear (List.tabulate (n, fn i => i))
    end

  (* ------------------------------------------------------------------ *)
  (* collide dispatch — manifolds *)

  fun aabbCollide (a : aabb) (b : aabb) : manifold option =
    let val a = normBox a and b = normBox b
    in if not (aabbOverlap a b) then NONE
       else
         let
           (* overlap on each axis *)
           val ox = fmin (V.x (#max a), V.x (#max b)) - fmax (V.x (#min a), V.x (#min b))
           val oy = fmin (V.y (#max a), V.y (#max b)) - fmax (V.y (#min a), V.y (#min b))
           val ca = V.scale (0.5, V.add (#min a, #max a))
           val cb = V.scale (0.5, V.add (#min b, #max b))
           val contact = V.scale (0.5, V.add (ca, cb))
         in
           if ox <= oy then
             let val sign = if V.x cb >= V.x ca then 1.0 else ~1.0
             in SOME { normal = V.v (sign, 0.0), depth = ox, contacts = [contact] } end
           else
             let val sign = if V.y cb >= V.y ca then 1.0 else ~1.0
             in SOME { normal = V.v (0.0, sign), depth = oy, contacts = [contact] } end
         end
    end

  fun circleCollide (c1 : circle) (c2 : circle) : manifold option =
    let val d = V.sub (#center c2, #center c1)
        val dist = V.length d
        val rs = Real.abs (#radius c1) + Real.abs (#radius c2)
    in if dist > rs then NONE
       else
         let val normal =
                   if dist < 1.0E~12 then V.v (1.0, 0.0)
                   else V.scale (1.0 / dist, d)
             val depth = rs - dist
             val contact = V.add (#center c1, V.scale (Real.abs (#radius c1), normal))
         in SOME { normal = normal, depth = depth, contacts = [contact] } end
    end

  fun polyCollide (ptsA, ptsB) : manifold option =
    case satPolys (ptsA, ptsB) of
      NONE => NONE
    | SOME (normal, depth) =>
        let
          fun centroid pts =
            let val n = Real.fromInt (List.length pts)
                val sx = List.foldl (fn (p, s) => s + V.x p) 0.0 pts
                val sy = List.foldl (fn (p, s) => s + V.y p) 0.0 pts
            in V.v (sx / n, sy / n) end
          val contact = V.scale (0.5, V.add (centroid ptsA, centroid ptsB))
        in SOME { normal = normal, depth = depth, contacts = [contact] } end

  fun collide s1 s2 =
    case (s1, s2) of
      (AABB a, AABB b) => aabbCollide a b
    | (Circle a, Circle b) => circleCollide a b
    | _ =>
        (case (shapeToPoly s1, shapeToPoly s2) of
           (SOME pa, SOME pb) => polyCollide (pa, pb)
         | _ => if overlap s1 s2
                then SOME { normal = V.v (0.0, 0.0), depth = 0.0, contacts = [] }
                else NONE)

  (* ------------------------------------------------------------------ *)
  (* Grid — pure spatial hash over fixed-size square cells. *)

  structure Grid =
  struct
    type cell = int * int
    type t = { cell : real,
               items : (int * shape) list,           (* id -> shape *)
               cells : (cell * int list) list }       (* cell -> ids *)

    fun make c = { cell = if c <= 0.0 then 1.0 else c, items = [], cells = [] }

    fun cellOf size x = Real.floor (x / size)

    (* All cells covered by a shape's AABB. *)
    fun cellsOfShape size sh =
      let val b = aabbOf sh
          val cx0 = cellOf size (V.x (#min b))
          val cy0 = cellOf size (V.y (#min b))
          val cx1 = cellOf size (V.x (#max b))
          val cy1 = cellOf size (V.y (#max b))
          val xs = List.tabulate (cx1 - cx0 + 1, fn i => cx0 + i)
          val ys = List.tabulate (cy1 - cy0 + 1, fn j => cy0 + j)
      in List.concat (List.map (fn ix => List.map (fn iy => (ix, iy)) ys) xs)
      end

    fun addToCell (cells, key, id) =
      let fun go [] = [(key, [id])]
            | go ((k, ids) :: rest) =
                if k = key then (k, id :: ids) :: rest
                else (k, ids) :: go rest
      in go cells end

    fun removeId (cells, id) =
      List.mapPartial
        (fn (k, ids) =>
           case List.filter (fn i => i <> id) ids of
             [] => NONE
           | kept => SOME (k, kept))
        cells

    fun remove (g : t) id =
      { cell = #cell g,
        items = List.filter (fn (i, _) => i <> id) (#items g),
        cells = removeId (#cells g, id) }

    fun insert (g : t) id sh =
      let val g = remove g id   (* re-insert semantics *)
          val keys = cellsOfShape (#cell g) sh
          val cells' = List.foldl (fn (k, cs) => addToCell (cs, k, id)) (#cells g) keys
      in { cell = #cell g, items = (id, sh) :: #items g, cells = cells' } end

    fun clear (g : t) = { cell = #cell g, items = [], cells = [] }

    fun query (g : t) sh =
      let val keys = cellsOfShape (#cell g) sh
          fun idsIn key =
            case List.find (fn (k, _) => k = key) (#cells g) of
              SOME (_, ids) => ids
            | NONE => []
          val raw = List.concat (List.map idsIn keys)
          (* dedup preserving order *)
          fun dedup ([], _) = []
            | dedup (x :: xs, seen) =
                if List.exists (fn s => s = x) seen then dedup (xs, seen)
                else x :: dedup (xs, x :: seen)
      in dedup (raw, []) end
  end
end
