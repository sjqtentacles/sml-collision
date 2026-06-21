(* glm.sml

   Implementation of GLM: pure linear algebra and transforms.

   Storage: vectors are tuples; matrices are tuples of their column-major
   elements (m00 is row 0 / col 0). Using plain `real` tuples keeps behaviour
   byte-identical across MLton and Poly/ML. *)

structure Glm :> GLM =
struct
  val pi = 3.14159265358979323846

  fun radians d = d * pi / 180.0
  fun degrees r = r * 180.0 / pi
  fun clamp (x, lo, hi) : real = if x < lo then lo else if x > hi then hi else x
  fun lerp (a, b, t) : real = a + (b - a) * t

  fun rstr r = Real.fmt (StringCvt.GEN (SOME 6)) r

  structure Vec2 =
  struct
    type t = real * real
    fun v p = p
    fun x (a, _) = a
    fun y (_, b) = b
    val zero = (0.0, 0.0)
    fun toList (a, b) = [a, b]

    fun add ((a,b),(c,d)) : t = (a+c, b+d)
    fun sub ((a,b),(c,d)) : t = (a-c, b-d)
    fun mulc ((a,b),(c,d)) : t = (a*c, b*d)
    fun scale (s,(a,b)) : t = (s*a, s*b)
    fun neg (a,b) : t = (~a, ~b)

    fun dot ((a,b),(c,d)) : real = a*c + b*d
    fun lengthSq u = dot (u, u)
    fun length u = Math.sqrt (lengthSq u)
    fun dist (u, w) = length (sub (u, w))
    fun normalize u =
      let val l = length u
      in if Real.== (l, 0.0) then zero else scale (1.0 / l, u) end
    fun lerp ((a,b),(c,d),t) : t = (a+(c-a)*t, b+(d-b)*t)

    fun equal ((a,b),(c,d)) = Real.== (a,c) andalso Real.== (b,d)
    fun approx eps ((a,b),(c,d)) =
      Real.abs (a-c) <= eps andalso Real.abs (b-d) <= eps
    fun toString (a,b) = "(" ^ rstr a ^ ", " ^ rstr b ^ ")"
  end

  structure Vec3 =
  struct
    type t = real * real * real
    fun v p = p
    fun x (a,_,_) = a
    fun y (_,b,_) = b
    fun z (_,_,c) = c
    val zero = (0.0, 0.0, 0.0)
    fun toList (a,b,c) = [a,b,c]

    fun add ((a,b,c),(d,e,f)) : t = (a+d, b+e, c+f)
    fun sub ((a,b,c),(d,e,f)) : t = (a-d, b-e, c-f)
    fun mulc ((a,b,c),(d,e,f)) : t = (a*d, b*e, c*f)
    fun scale (s,(a,b,c)) : t = (s*a, s*b, s*c)
    fun neg (a,b,c) : t = (~a, ~b, ~c)

    fun dot ((a,b,c),(d,e,f)) : real = a*d + b*e + c*f
    fun cross ((a,b,c),(d,e,f)) : t = (b*f - c*e, c*d - a*f, a*e - b*d)
    fun lengthSq u = dot (u, u)
    fun length u = Math.sqrt (lengthSq u)
    fun dist (u, w) = length (sub (u, w))
    fun normalize u =
      let val l = length u
      in if Real.== (l, 0.0) then zero else scale (1.0 / l, u) end
    fun lerp ((a,b,c),(d,e,f),t) : t = (a+(d-a)*t, b+(e-b)*t, c+(f-c)*t)

    fun equal ((a,b,c),(d,e,f)) =
      Real.== (a,d) andalso Real.== (b,e) andalso Real.== (c,f)
    fun approx eps ((a,b,c),(d,e,f)) =
      Real.abs (a-d) <= eps andalso Real.abs (b-e) <= eps
      andalso Real.abs (c-f) <= eps
    fun toString (a,b,c) = "(" ^ rstr a ^ ", " ^ rstr b ^ ", " ^ rstr c ^ ")"
  end

  structure Vec4 =
  struct
    type t = real * real * real * real
    fun v p = p
    fun x (a,_,_,_) = a
    fun y (_,b,_,_) = b
    fun z (_,_,c,_) = c
    fun w (_,_,_,d) = d
    val zero = (0.0, 0.0, 0.0, 0.0)
    fun toList (a,b,c,d) = [a,b,c,d]

    fun add ((a,b,c,d),(e,f,g,h)) : t = (a+e, b+f, c+g, d+h)
    fun sub ((a,b,c,d),(e,f,g,h)) : t = (a-e, b-f, c-g, d-h)
    fun mulc ((a,b,c,d),(e,f,g,h)) : t = (a*e, b*f, c*g, d*h)
    fun scale (s,(a,b,c,d)) : t = (s*a, s*b, s*c, s*d)
    fun neg (a,b,c,d) : t = (~a, ~b, ~c, ~d)

    fun dot ((a,b,c,d),(e,f,g,h)) : real = a*e + b*f + c*g + d*h
    fun lengthSq u = dot (u, u)
    fun length u = Math.sqrt (lengthSq u)
    fun dist (u, ww) = length (sub (u, ww))
    fun normalize u =
      let val l = length u
      in if Real.== (l, 0.0) then zero else scale (1.0 / l, u) end
    fun lerp ((a,b,c,d),(e,f,g,h),t) : t =
      (a+(e-a)*t, b+(f-b)*t, c+(g-c)*t, d+(h-d)*t)

    fun equal ((a,b,c,d),(e,f,g,h)) =
      Real.== (a,e) andalso Real.== (b,f) andalso Real.== (c,g)
      andalso Real.== (d,h)
    fun approx eps ((a,b,c,d),(e,f,g,h)) =
      Real.abs (a-e) <= eps andalso Real.abs (b-f) <= eps
      andalso Real.abs (c-g) <= eps andalso Real.abs (d-h) <= eps
    fun toString (a,b,c,d) =
      "(" ^ rstr a ^ ", " ^ rstr b ^ ", " ^ rstr c ^ ", " ^ rstr d ^ ")"
  end

  fun matToString (label, els) =
    label ^ "[" ^ String.concatWith ", " (List.map rstr els) ^ "]"

  structure Mat3 =
  struct
    (* Column-major: (m00,m10,m20, m01,m11,m21, m02,m12,m22) where mRC is
       row R, column C. *)
    type t = real * real * real * real * real * real * real * real * real

    val id = (1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0)

    fun fromCols ((a,b,c),(d,e,f),(g,h,i)) = (a,b,c,d,e,f,g,h,i)
    fun fromRows ((a,b,c),(d,e,f),(g,h,i)) =
      (* rows -> column-major *)
      (a,d,g, b,e,h, c,f,i)
    fun fromList [a,b,c,d,e,f,g,h,i] = (a,b,c,d,e,f,g,h,i)
      | fromList _ = raise Fail "Mat3.fromList: need 9 elements"
    fun toList (a,b,c,d,e,f,g,h,i) = [a,b,c,d,e,f,g,h,i]

    fun add ((a,b,c,d,e,f,g,h,i),(a',b',c',d',e',f',g',h',i')) : t =
      (a+a',b+b',c+c',d+d',e+e',f+f',g+g',h+h',i+i')
    fun sub ((a,b,c,d,e,f,g,h,i),(a',b',c',d',e',f',g',h',i')) : t =
      (a-a',b-b',c-c',d-d',e-e',f-f',g-g',h-h',i-i')
    fun scale (s,(a,b,c,d,e,f,g,h,i)) : t =
      (s*a,s*b,s*c,s*d,s*e,s*f,s*g,s*h,s*i)

    (* element access by (row, col), column-major storage *)
    (* m = (m00,m10,m20, m01,m11,m21, m02,m12,m22) *)
    (* These matrix routines route their many real intermediates through a
     * mutable array rather than binding them all as simultaneously-live locals.
     * That keeps floating-point register pressure low enough to dodge Poly/ML's
     * native-codegen bugs ("asFPReg"/"getAllocatedGenReg raised while
     * compiling") on large straight-line real arithmetic. Behaviour is
     * identical to the obvious tuple version. *)
    fun mul ((a00,a10,a20, a01,a11,a21, a02,a12,a22) : t,
             (b00,b10,b20, b01,b11,b21, b02,b12,b22) : t) : t =
      let
        val r = Array.array (9, 0.0)
        (* result row i, col j = sum_k A[i,k] B[k,j]; stored column-major. *)
        val () = Array.update (r, 0, a00*b00 + a01*b10 + a02*b20)  (* r00 *)
        val () = Array.update (r, 1, a10*b00 + a11*b10 + a12*b20)  (* r10 *)
        val () = Array.update (r, 2, a20*b00 + a21*b10 + a22*b20)  (* r20 *)
        val () = Array.update (r, 3, a00*b01 + a01*b11 + a02*b21)  (* r01 *)
        val () = Array.update (r, 4, a10*b01 + a11*b11 + a12*b21)  (* r11 *)
        val () = Array.update (r, 5, a20*b01 + a21*b11 + a22*b21)  (* r21 *)
        val () = Array.update (r, 6, a00*b02 + a01*b12 + a02*b22)  (* r02 *)
        val () = Array.update (r, 7, a10*b02 + a11*b12 + a12*b22)  (* r12 *)
        val () = Array.update (r, 8, a20*b02 + a21*b12 + a22*b22)  (* r22 *)
        fun g i = Array.sub (r, i)
      in
        (g 0, g 1, g 2, g 3, g 4, g 5, g 6, g 7, g 8)
      end

    fun transpose ((m00,m10,m20, m01,m11,m21, m02,m12,m22) : t) : t =
      (m00,m01,m02, m10,m11,m12, m20,m21,m22)

    fun det ((m00,m10,m20, m01,m11,m21, m02,m12,m22) : t) : real =
      m00*(m11*m22 - m12*m21)
      - m01*(m10*m22 - m12*m20)
      + m02*(m10*m21 - m11*m20)

    fun inverse (m as (m00,m10,m20, m01,m11,m21, m02,m12,m22) : t) : t option =
      let
        val d = det m
      in
        if Real.== (d, 0.0) then NONE
        else
          let
            val inv = 1.0 / d
            (* cofactors into an array, then inverse[r,c] = inv * cofactor[c,r];
               result stored column-major. *)
            val c = Array.array (9, 0.0)
            val () = Array.update (c, 0,  (m11*m22 - m12*m21))  (* c00 *)
            val () = Array.update (c, 1, ~(m10*m22 - m12*m20))  (* c01 *)
            val () = Array.update (c, 2,  (m10*m21 - m11*m20))  (* c02 *)
            val () = Array.update (c, 3, ~(m01*m22 - m02*m21))  (* c10 *)
            val () = Array.update (c, 4,  (m00*m22 - m02*m20))  (* c11 *)
            val () = Array.update (c, 5, ~(m00*m21 - m01*m20))  (* c12 *)
            val () = Array.update (c, 6,  (m01*m12 - m02*m11))  (* c20 *)
            val () = Array.update (c, 7, ~(m00*m12 - m02*m10))  (* c21 *)
            val () = Array.update (c, 8,  (m00*m11 - m01*m10))  (* c22 *)
            fun co i = inv * Array.sub (c, i)
          in
            (* (i00,i10,i20, i01,i11,i21, i02,i12,i22)
               where i[r,c] = inv*cofactor[c,r] *)
            SOME (co 0, co 1, co 2, co 3, co 4, co 5, co 6, co 7, co 8)
          end
      end

    fun mulV ((m00,m10,m20, m01,m11,m21, m02,m12,m22) : t,(x,y,z)) : Vec3.t =
      (m00*x + m01*y + m02*z,
       m10*x + m11*y + m12*z,
       m20*x + m21*y + m22*z)

    fun equal (a, b) =
      ListPair.allEq (fn (p,q) => Real.== (p,q)) (toList a, toList b)
    fun approx eps (a, b) =
      ListPair.all (fn (p,q) => Real.abs (p-q) <= eps) (toList a, toList b)
    fun toString m = matToString ("Mat3", toList m)
  end

  structure Mat4 =
  struct
    (* Column-major, 16 elements in a vector. Index (row r, col c) = c*4 + r. *)
    type t = real vector

    fun idx (m : t, r, c) : real = Vector.sub (m, c*4 + r)

    val id = Vector.fromList
      [1.0,0.0,0.0,0.0,
       0.0,1.0,0.0,0.0,
       0.0,0.0,1.0,0.0,
       0.0,0.0,0.0,1.0]

    fun fromList xs =
      if List.length xs = 16 then Vector.fromList xs
      else raise Fail "Mat4.fromList: need 16 elements"
    fun toList m = Vector.foldr (op ::) [] m

    fun fromCols ((a,b,c,d),(e,f,g,h),(i,j,k,l),(m,n,p,q)) =
      Vector.fromList [a,b,c,d, e,f,g,h, i,j,k,l, m,n,p,q]
    fun fromRows ((a,b,c,d),(e,f,g,h),(i,j,k,l),(m,n,p,q)) =
      (* rows -> column major *)
      Vector.fromList [a,e,i,m, b,f,j,n, c,g,k,p, d,h,l,q]

    fun binop (f : real * real -> real) (a, b) =
      Vector.tabulate (16, fn k => f (Vector.sub (a,k), Vector.sub (b,k)))
    fun add p = binop (op +) p
    fun sub p = binop (op -) p
    fun scale (s : real, m) = Vector.map (fn e => s*e) m

    fun mul (a, b) =
      Vector.tabulate (16, fn k =>
        let val c = k div 4   (* column of result *)
            val r = k mod 4   (* row of result *)
        in idx (a,r,0)*idx (b,0,c) + idx (a,r,1)*idx (b,1,c)
           + idx (a,r,2)*idx (b,2,c) + idx (a,r,3)*idx (b,3,c)
        end)

    fun transpose m =
      Vector.tabulate (16, fn k =>
        let val c = k div 4 val r = k mod 4 in idx (m, c, r) end)

    fun mulV (m, (x,y,z,w)) =
      (idx(m,0,0)*x + idx(m,0,1)*y + idx(m,0,2)*z + idx(m,0,3)*w,
       idx(m,1,0)*x + idx(m,1,1)*y + idx(m,1,2)*z + idx(m,1,3)*w,
       idx(m,2,0)*x + idx(m,2,1)*y + idx(m,2,2)*z + idx(m,2,3)*w,
       idx(m,3,0)*x + idx(m,3,1)*y + idx(m,3,2)*z + idx(m,3,3)*w)

    fun transformPoint (m, (x,y,z)) =
      let val (rx,ry,rz,rw) = mulV (m, (x,y,z,1.0))
      in if Real.== (rw, 0.0) then (rx,ry,rz)
         else (rx/rw, ry/rw, rz/rw) end
    fun transformDir (m, (x,y,z)) =
      let val (rx,ry,rz,_) = mulV (m, (x,y,z,0.0)) in (rx,ry,rz) end

    (* Laplace expansion for determinant and inverse via adjugate. *)
    fun det m =
      let
        fun e (r,c) = idx (m,r,c)
        val s0 = e(0,0)*e(1,1) - e(1,0)*e(0,1)
        val s1 = e(0,0)*e(1,2) - e(1,0)*e(0,2)
        val s2 = e(0,0)*e(1,3) - e(1,0)*e(0,3)
        val s3 = e(0,1)*e(1,2) - e(1,1)*e(0,2)
        val s4 = e(0,1)*e(1,3) - e(1,1)*e(0,3)
        val s5 = e(0,2)*e(1,3) - e(1,2)*e(0,3)
        val c5 = e(2,2)*e(3,3) - e(3,2)*e(2,3)
        val c4 = e(2,1)*e(3,3) - e(3,1)*e(2,3)
        val c3 = e(2,1)*e(3,2) - e(3,1)*e(2,2)
        val c2 = e(2,0)*e(3,3) - e(3,0)*e(2,3)
        val c1 = e(2,0)*e(3,2) - e(3,0)*e(2,2)
        val c0 = e(2,0)*e(3,1) - e(3,0)*e(2,1)
      in
        s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0
      end

    fun inverse m =
      let
        fun e (r,c) = idx (m,r,c)
        val s0 = e(0,0)*e(1,1) - e(1,0)*e(0,1)
        val s1 = e(0,0)*e(1,2) - e(1,0)*e(0,2)
        val s2 = e(0,0)*e(1,3) - e(1,0)*e(0,3)
        val s3 = e(0,1)*e(1,2) - e(1,1)*e(0,2)
        val s4 = e(0,1)*e(1,3) - e(1,1)*e(0,3)
        val s5 = e(0,2)*e(1,3) - e(1,2)*e(0,3)
        val c5 = e(2,2)*e(3,3) - e(3,2)*e(2,3)
        val c4 = e(2,1)*e(3,3) - e(3,1)*e(2,3)
        val c3 = e(2,1)*e(3,2) - e(3,1)*e(2,2)
        val c2 = e(2,0)*e(3,3) - e(3,0)*e(2,3)
        val c1 = e(2,0)*e(3,2) - e(3,0)*e(2,2)
        val c0 = e(2,0)*e(3,1) - e(3,0)*e(2,1)
        val d = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0
      in
        if Real.== (d, 0.0) then NONE
        else
          let
            val invd = 1.0 / d
            fun ne (r,c) = ~(e (r,c))
            (* Compute the 16 entries one at a time straight into a mutable
               array (row-major here), so they are never all simultaneously
               live. This keeps Poly/ML's native codegen from exhausting FP
               registers ("asFPReg raised while compiling"). *)
            val b = Array.array (16, 0.0)
            fun set (i, v) = Array.update (b, i, v * invd)
            val () = set (0,   e(1,1)*c5 - e(1,2)*c4 + e(1,3)*c3)
            val () = set (1,  ne(0,1)*c5 + e(0,2)*c4 - e(0,3)*c3)
            val () = set (2,   e(3,1)*s5 - e(3,2)*s4 + e(3,3)*s3)
            val () = set (3,  ne(2,1)*s5 + e(2,2)*s4 - e(2,3)*s3)
            val () = set (4,  ne(1,0)*c5 + e(1,2)*c2 - e(1,3)*c1)
            val () = set (5,   e(0,0)*c5 - e(0,2)*c2 + e(0,3)*c1)
            val () = set (6,  ne(3,0)*s5 + e(3,2)*s2 - e(3,3)*s1)
            val () = set (7,   e(2,0)*s5 - e(2,2)*s2 + e(2,3)*s1)
            val () = set (8,   e(1,0)*c4 - e(1,1)*c2 + e(1,3)*c0)
            val () = set (9,  ne(0,0)*c4 + e(0,1)*c2 - e(0,3)*c0)
            val () = set (10,  e(3,0)*s4 - e(3,1)*s2 + e(3,3)*s0)
            val () = set (11, ne(2,0)*s4 + e(2,1)*s2 - e(2,3)*s0)
            val () = set (12, ne(1,0)*c3 + e(1,1)*c1 - e(1,2)*c0)
            val () = set (13,  e(0,0)*c3 - e(0,1)*c1 + e(0,2)*c0)
            val () = set (14, ne(3,0)*s3 + e(3,1)*s1 - e(3,2)*s0)
            val () = set (15,  e(2,0)*s3 - e(2,1)*s1 + e(2,2)*s0)
            fun g (r, c) = Array.sub (b, r*4 + c)
          in
            SOME (fromRows ((g(0,0),g(0,1),g(0,2),g(0,3)),
                            (g(1,0),g(1,1),g(1,2),g(1,3)),
                            (g(2,0),g(2,1),g(2,2),g(2,3)),
                            (g(3,0),g(3,1),g(3,2),g(3,3))))
          end
      end

    fun translate (x,y,z) =
      fromRows ((1.0,0.0,0.0,x),(0.0,1.0,0.0,y),
                (0.0,0.0,1.0,z),(0.0,0.0,0.0,1.0))
    fun scaleM (x,y,z) =
      fromRows ((x,0.0,0.0,0.0),(0.0,y,0.0,0.0),
                (0.0,0.0,z,0.0),(0.0,0.0,0.0,1.0))

    fun rotate (angle, axis) =
      let
        val (ax,ay,az) = Vec3.normalize axis
        val c = Math.cos angle
        val s = Math.sin angle
        val t = 1.0 - c
      in
        fromRows
          (( t*ax*ax + c,      t*ax*ay - s*az,  t*ax*az + s*ay,  0.0),
           ( t*ax*ay + s*az,   t*ay*ay + c,     t*ay*az - s*ax,  0.0),
           ( t*ax*az - s*ay,   t*ay*az + s*ax,  t*az*az + c,     0.0),
           ( 0.0,              0.0,             0.0,             1.0))
      end
    fun rotateX a = rotate (a, (1.0,0.0,0.0))
    fun rotateY a = rotate (a, (0.0,1.0,0.0))
    fun rotateZ a = rotate (a, (0.0,0.0,1.0))

    (* Build a Mat4 from a function giving each column-major element by index.
       Routing every builder through this single, separately-compiled function
       (rather than inlining a 16-wide tuple + Vector.fromList at each call
       site) keeps Poly/ML's native codegen from exhausting its registers
       ("asGenReg"/"asFPReg raised while compiling") on older x86-64 builds. *)
    fun build (f : int -> real) : t = Vector.tabulate (16, f)

    fun perspective {fovy, aspect, near, far} =
      let
        val ft = 1.0 / Math.tan (fovy / 2.0)
        val nf = 1.0 / (near - far)
        val a = Array.array (16, 0.0)
        val () = Array.update (a, 0, ft / aspect)
        val () = Array.update (a, 5, ft)
        val () = Array.update (a, 10, (far + near) * nf)
        val () = Array.update (a, 11, ~1.0)
        val () = Array.update (a, 14, 2.0 * far * near * nf)
      in
        build (fn k => Array.sub (a, k))
      end

    fun ortho {left, right, bottom, top, near, far} =
      let
        val rl = 1.0 / (right - left)
        val tb = 1.0 / (top - bottom)
        val fn_ = 1.0 / (far - near)
        val a = Array.array (16, 0.0)
        val () = Array.update (a, 0, 2.0 * rl)
        val () = Array.update (a, 5, 2.0 * tb)
        val () = Array.update (a, 10, ~2.0 * fn_)
        val () = Array.update (a, 12, ~(right + left) * rl)
        val () = Array.update (a, 13, ~(top + bottom) * tb)
        val () = Array.update (a, 14, ~(far + near) * fn_)
        val () = Array.update (a, 15, 1.0)
      in
        build (fn k => Array.sub (a, k))
      end

    fun lookAt {eye, center, up} =
      let
        val f = Vec3.normalize (Vec3.sub (center, eye))   (* forward *)
        val s = Vec3.normalize (Vec3.cross (f, up))       (* right *)
        val u = Vec3.cross (s, f)                         (* true up *)
        val (sx,sy,sz) = s
        val (ux,uy,uz) = u
        val (fx,fy,fz) = f
        (* Build the 16 entries one at a time through a mutable array so they
           are never all simultaneously live; this avoids Poly/ML's native
           codegen FP-register exhaustion ("asFPReg raised while compiling"). *)
        val b = Array.array (16, 0.0)
        fun set (i, v) = Array.update (b, i, v)
        val () = set (0, sx)  val () = set (1, sy)  val () = set (2, sz)
        val () = set (3, ~(Vec3.dot (s, eye)))
        val () = set (4, ux)  val () = set (5, uy)  val () = set (6, uz)
        val () = set (7, ~(Vec3.dot (u, eye)))
        val () = set (8, ~fx) val () = set (9, ~fy) val () = set (10, ~fz)
        val () = set (11, Vec3.dot (f, eye))
        val () = set (12, 0.0) val () = set (13, 0.0)
        val () = set (14, 0.0) val () = set (15, 1.0)
        fun g (r, c) = Array.sub (b, r*4 + c)
      in
        fromRows ((g(0,0),g(0,1),g(0,2),g(0,3)),
                  (g(1,0),g(1,1),g(1,2),g(1,3)),
                  (g(2,0),g(2,1),g(2,2),g(2,3)),
                  (g(3,0),g(3,1),g(3,2),g(3,3)))
      end

    fun equal (a, b) =
      ListPair.allEq (fn (p,q) => Real.== (p,q)) (toList a, toList b)
    fun approx eps (a, b) =
      ListPair.all (fn (p,q) => Real.abs (p-q) <= eps) (toList a, toList b)
    fun toString m = matToString ("Mat4", toList m)
  end

  structure Quat =
  struct
    (* (w, x, y, z) *)
    type t = real * real * real * real
    fun quat p = p
    val id = (1.0, 0.0, 0.0, 0.0)
    fun w (a,_,_,_) = a
    fun x (_,b,_,_) = b
    fun y (_,_,c,_) = c
    fun z (_,_,_,d) = d
    fun toList (a,b,c,d) = [a,b,c,d]

    fun fromAxisAngle (axis, angle) =
      let
        val (ax,ay,az) = Vec3.normalize axis
        val half = angle / 2.0
        val s = Math.sin half
      in
        (Math.cos half, ax*s, ay*s, az*s)
      end

    fun add ((a,b,c,d),(e,f,g,h)) : t = (a+e, b+f, c+g, d+h)
    fun scale (s,(a,b,c,d)) : t = (s*a, s*b, s*c, s*d)

    (* Hamilton product *)
    fun mul ((w1,x1,y1,z1),(w2,x2,y2,z2)) : t =
      (w1*w2 - x1*x2 - y1*y2 - z1*z2,
       w1*x2 + x1*w2 + y1*z2 - z1*y2,
       w1*y2 - x1*z2 + y1*w2 + z1*x2,
       w1*z2 + x1*y2 - y1*x2 + z1*w2)

    fun conj (a,b,c,d) : t = (a, ~b, ~c, ~d)
    fun dot ((a,b,c,d),(e,f,g,h)) : real = a*e + b*f + c*g + d*h
    fun length q = Math.sqrt (dot (q, q))
    fun normalize q =
      let val l = length q
      in if Real.== (l, 0.0) then id else scale (1.0 / l, q) end

    fun rotateV (q, (vx,vy,vz)) =
      let
        val (_, rx, ry, rz) = mul (mul (q, (0.0, vx, vy, vz)), conj q)
      in
        (rx, ry, rz)
      end

    fun sub (a, b) = add (a, scale (~1.0, b))

    (* Factored out so the trig/coefficient computation compiles as its own
       code unit; folding it inline pushes Poly/ML's native codegen over its
       FP-register budget ("asFPReg raised while compiling"). *)
    fun slerpCoeffs (d, t) =
      let
        val theta0 = Math.acos (clamp (d, ~1.0, 1.0))
        val theta = theta0 * t
        val sinTheta0 = Math.sin theta0
        val s2 = Math.sin theta / sinTheta0
        val s1 = Math.cos theta - d * s2
      in
        (s1, s2)
      end

    (* The blend is its own top-level function (not inlined into slerp) and
       writes through a mutable array, so neither slerp nor the blend holds the
       full set of quaternion component reals live at once. This avoids
       Poly/ML's native-codegen register exhaustion
       ("asGenReg"/"asFPReg raised while compiling") on older x86-64 builds. *)
    fun blend (c1, c2, (w1,x1,y1,z1), (w2,x2,y2,z2)) : t =
      let
        val r = Array.array (4, 0.0)
        val () = Array.update (r, 0, c1*w1 + c2*w2)
        val () = Array.update (r, 1, c1*x1 + c2*x2)
        val () = Array.update (r, 2, c1*y1 + c2*y2)
        val () = Array.update (r, 3, c1*z1 + c2*z2)
      in
        (Array.sub (r,0), Array.sub (r,1), Array.sub (r,2), Array.sub (r,3))
      end

    fun slerp (q1, q2, t) =
      let
        val d0 = dot (q1, q2)
        (* take shorter arc *)
        val (q2, d) = if d0 < 0.0 then (scale (~1.0, q2), ~d0) else (q2, d0)
        val near = d > 0.9995
        val (c1, c2) =
          if near then (1.0 - t, t)  (* nearly parallel: linear *)
          else slerpCoeffs (d, t)
        val res = blend (c1, c2, q1, q2)
      in
        (* The linear branch must be renormalised; the slerp branch is already
           unit-length up to rounding but normalising is harmless. *)
        if near then normalize res else res
      end

    fun toMat4 q =
      let
        val (qw,qx,qy,qz) = normalize q
        (* Build the rotation matrix one entry at a time through a mutable
           array so the many products are never all simultaneously live; this
           avoids Poly/ML's native-codegen register exhaustion
           ("asFPReg"/"getAllocatedGenReg raised while compiling"). *)
        val xx = qx*qx val yy = qy*qy val zz = qz*qz
        val xy = qx*qy val xz = qx*qz val yz = qy*qz
        val wx = qw*qx val wy = qw*qy val wz = qw*qz
        val b = Array.array (16, 0.0)
        fun set (i, v) = Array.update (b, i, v)
        val () = set (0, 1.0 - 2.0*(yy+zz))
        val () = set (1, 2.0*(xy - wz))
        val () = set (2, 2.0*(xz + wy))
        val () = set (4, 2.0*(xy + wz))
        val () = set (5, 1.0 - 2.0*(xx+zz))
        val () = set (6, 2.0*(yz - wx))
        val () = set (8, 2.0*(xz - wy))
        val () = set (9, 2.0*(yz + wx))
        val () = set (10, 1.0 - 2.0*(xx+yy))
        val () = set (15, 1.0)
        fun g (r, c) = Array.sub (b, r*4 + c)
      in
        Mat4.fromRows ((g(0,0),g(0,1),g(0,2),g(0,3)),
                       (g(1,0),g(1,1),g(1,2),g(1,3)),
                       (g(2,0),g(2,1),g(2,2),g(2,3)),
                       (g(3,0),g(3,1),g(3,2),g(3,3)))
      end

    fun fromMat3 m =
      let
        val e = Mat3.toList m
        (* column-major (m00,m10,m20, m01,m11,m21, m02,m12,m22) *)
        fun el (r, c) = List.nth (e, c*3 + r)
        val tr = el(0,0) + el(1,1) + el(2,2)
      in
        if tr > 0.0 then
          let
            val s = Math.sqrt (tr + 1.0) * 2.0   (* s = 4*qw *)
          in
            (0.25*s,
             (el(2,1) - el(1,2)) / s,
             (el(0,2) - el(2,0)) / s,
             (el(1,0) - el(0,1)) / s)
          end
        else if el(0,0) > el(1,1) andalso el(0,0) > el(2,2) then
          let val s = Math.sqrt (1.0 + el(0,0) - el(1,1) - el(2,2)) * 2.0
          in ((el(2,1) - el(1,2)) / s, 0.25*s,
              (el(0,1) + el(1,0)) / s, (el(0,2) + el(2,0)) / s) end
        else if el(1,1) > el(2,2) then
          let val s = Math.sqrt (1.0 + el(1,1) - el(0,0) - el(2,2)) * 2.0
          in ((el(0,2) - el(2,0)) / s, (el(0,1) + el(1,0)) / s,
              0.25*s, (el(1,2) + el(2,1)) / s) end
        else
          let val s = Math.sqrt (1.0 + el(2,2) - el(0,0) - el(1,1)) * 2.0
          in ((el(1,0) - el(0,1)) / s, (el(0,2) + el(2,0)) / s,
              (el(1,2) + el(2,1)) / s, 0.25*s) end
      end

    fun equal ((a,b,c,d),(e,f,g,h)) =
      Real.== (a,e) andalso Real.== (b,f) andalso Real.== (c,g)
      andalso Real.== (d,h)
    fun approx eps ((a,b,c,d),(e,f,g,h)) =
      Real.abs (a-e) <= eps andalso Real.abs (b-f) <= eps
      andalso Real.abs (c-g) <= eps andalso Real.abs (d-h) <= eps
    fun toString (a,b,c,d) =
      "Quat(" ^ rstr a ^ ", " ^ rstr b ^ ", " ^ rstr c ^ ", " ^ rstr d ^ ")"
  end
end
