(* test.sml — tests for sml-collision. *)
structure Tests =
struct
  open Collision
  structure V = Collision.Glm.Vec2

  (* local real-approx helper *)
  val eps = 1.0E~6
  fun checkReal name (expected, actual) =
    Harness.check (name ^ " (" ^ Real.toString expected ^ " ~ " ^ Real.toString actual ^ ")")
                  (Real.abs (expected - actual) < eps)

  fun v (x, y) = V.v (x, y)

  (* a unit AABB with corner at (ox,oy) and size 1x1 *)
  fun box (ox, oy) = AABB { min = v (ox, oy), max = v (ox + 1.0, oy + 1.0) }
  fun circ (cx, cy, r) = Circle { center = v (cx, cy), radius = r }

  (* unit square as a CCW polygon with lower-left corner at (ox,oy) *)
  fun sqPoly (ox, oy) =
    Poly [ v (ox, oy), v (ox + 1.0, oy), v (ox + 1.0, oy + 1.0), v (ox, oy + 1.0) ]

  fun run () =
    let
      val () = Harness.reset ()

      (* ---------------------------------------------------------------- *)
      val () = Harness.section "AABB vs AABB overlap"
      val () = Harness.checkBool "overlapping (offset 0.5)"
                 (true, overlap (box (0.0, 0.0)) (box (0.5, 0.5)))
      val () = Harness.checkBool "touching exactly (offset 1.0) -> closed true"
                 (true, overlap (box (0.0, 0.0)) (box (1.0, 0.0)))
      val () = Harness.checkBool "separated (offset 2.0)"
                 (false, overlap (box (0.0, 0.0)) (box (2.0, 0.0)))
      val () = Harness.checkBool "separated on y only"
                 (false, overlap (box (0.0, 0.0)) (box (0.5, 2.0)))

      (* ---------------------------------------------------------------- *)
      val () = Harness.section "AABB vs AABB collide"
      val () =
        case collide (box (0.0, 0.0)) (box (0.5, 0.0)) of
          NONE => Harness.check "offset 0.5 on x -> SOME manifold" false
        | SOME { normal, depth, contacts } =>
            (Harness.check "depth > 0" (depth > 0.0);
             checkReal "depth ~ 0.5" (0.5, depth);
             Harness.check "normal points along x"
               (Real.abs (V.x normal) > Real.abs (V.y normal));
             Harness.check "has a contact point" (not (List.null contacts)))
      val () = Harness.checkBool "separated boxes -> NONE"
                 (true, not (Option.isSome (collide (box (0.0,0.0)) (box (3.0,0.0)))))
      val () =
        case collide (box (0.0, 0.0)) (box (0.0, 0.4)) of
          NONE => Harness.check "y-overlap -> SOME manifold" false
        | SOME { normal, depth, ... } =>
            (checkReal "depth ~ 0.6 on y axis" (0.6, depth);
             Harness.check "normal points along y"
               (Real.abs (V.y normal) > Real.abs (V.x normal)))

      (* ---------------------------------------------------------------- *)
      val () = Harness.section "Circle vs Circle overlap"
      val () = Harness.checkBool "concentric"
                 (true, overlap (circ (0.0,0.0,1.0)) (circ (0.0,0.0,0.5)))
      val () = Harness.checkBool "exactly r1+r2 apart -> boundary true"
                 (true, overlap (circ (0.0,0.0,1.0)) (circ (3.0,0.0,2.0)))
      val () = Harness.checkBool "farther apart -> false"
                 (false, overlap (circ (0.0,0.0,1.0)) (circ (3.1,0.0,2.0)))
      val () =
        case collide (circ (0.0,0.0,1.0)) (circ (1.5,0.0,1.0)) of
          NONE => Harness.check "overlapping circles -> SOME" false
        | SOME { normal, depth, ... } =>
            (checkReal "depth ~ 0.5" (0.5, depth);
             checkReal "normal is unit length" (1.0, V.length normal);
             Harness.check "normal points along +x" (V.x normal > 0.0))

      (* ---------------------------------------------------------------- *)
      val () = Harness.section "AABB vs Circle"
      val () = Harness.checkBool "circle center inside box"
                 (true, overlap (box (0.0,0.0)) (circ (0.5,0.5,0.1)))
      val () = Harness.checkBool "circle far outside box"
                 (false, overlap (box (0.0,0.0)) (circ (5.0,5.0,0.5)))
      val () = Harness.checkBool "circle tangent to box edge -> boundary true"
                 (true, overlap (box (0.0,0.0)) (circ (1.5,0.5,0.5)))
      val () = Harness.checkBool "circle just past edge -> false"
                 (false, overlap (box (0.0,0.0)) (circ (1.6,0.5,0.5)))
      val () = Harness.checkBool "symmetric: circle then box"
                 (true, overlap (circ (0.5,0.5,0.1)) (box (0.0,0.0)))

      (* ---------------------------------------------------------------- *)
      val () = Harness.section "Poly SAT"
      val tri = Poly [ v (0.0, 0.0), v (2.0, 0.0), v (1.0, 1.7) ]
      val triFar = Poly [ v (10.0, 10.0), v (12.0, 10.0), v (11.0, 11.7) ]
      val () = Harness.checkBool "triangle overlaps itself" (true, overlap tri tri)
      val () = Harness.checkBool "triangle vs far translate -> false"
                 (false, overlap tri triFar)
      val () = Harness.checkBool "two unit squares overlapping (offset 0.5)"
                 (true, overlap (sqPoly (0.0,0.0)) (sqPoly (0.5,0.5)))
      val () = Harness.checkBool "two unit squares with clear gap -> false"
                 (false, overlap (sqPoly (0.0,0.0)) (sqPoly (3.0,0.0)))
      val () = Harness.checkBool "poly square vs AABB overlapping"
                 (true, overlap (sqPoly (0.0,0.0)) (box (0.5,0.5)))
      val () = Harness.checkBool "poly square vs AABB separated"
                 (false, overlap (sqPoly (0.0,0.0)) (box (5.0,5.0)))
      val () =
        case collide (sqPoly (0.0,0.0)) (sqPoly (0.5,0.0)) of
          NONE => Harness.check "overlapping polys -> SOME manifold" false
        | SOME { depth, ... } =>
            (Harness.check "poly collide depth > 0" (depth > 0.0);
             checkReal "poly collide depth ~ 0.5" (0.5, depth))

      (* ---------------------------------------------------------------- *)
      val () = Harness.section "Segment"
      val seg = Segment { a = v (0.5, 0.5), b = v (0.5, 0.5) }
      val () = Harness.checkBool "degenerate segment inside box"
                 (true, overlap seg (box (0.0,0.0)))
      val () = Harness.checkBool "degenerate segment outside box"
                 (false, overlap (Segment { a = v (5.0,5.0), b = v (5.0,5.0) }) (box (0.0,0.0)))
      val () = Harness.checkBool "point inside circle"
                 (true, overlap (Segment { a = v (0.0,0.0), b = v (0.0,0.0) }) (circ (0.0,0.0,1.0)))
      val () = Harness.checkBool "segment crossing box"
                 (true, overlap (Segment { a = v (~1.0,0.5), b = v (2.0,0.5) }) (box (0.0,0.0)))

      (* ---------------------------------------------------------------- *)
      val () = Harness.section "Grid spatial hash"
      val g0 = Grid.make 1.0
      val g1 = Grid.insert g0 1 (box (0.0, 0.0))
      val g2 = Grid.insert g1 2 (box (5.0, 5.0))
      val g3 = Grid.insert g2 3 (box (10.0, 10.0))
      (* query near id 1 only *)
      val near1 = Grid.query g3 (box (0.0, 0.0))
      val () = Harness.checkBool "query near box 1 returns id 1"
                 (true, List.exists (fn i => i = 1) near1)
      val () = Harness.checkBool "query near box 1 excludes id 3"
                 (false, List.exists (fn i => i = 3) near1)
      val near2 = Grid.query g3 (box (5.0, 5.0))
      val () = Harness.checkBool "query near box 2 returns id 2"
                 (true, List.exists (fn i => i = 2) near2)
      val () = Harness.checkBool "query near box 2 excludes id 1"
                 (false, List.exists (fn i => i = 1) near2)
      (* removal *)
      val g4 = Grid.remove g3 2
      val () = Harness.checkBool "after remove, id 2 gone"
                 (false, List.exists (fn i => i = 2) (Grid.query g4 (box (5.0,5.0))))
      (* clear *)
      val g5 = Grid.clear g3
      val () = Harness.checkBool "after clear, query empty"
                 (true, List.null (Grid.query g5 (box (0.0,0.0))))
    in
      Harness.run ()
    end
end
