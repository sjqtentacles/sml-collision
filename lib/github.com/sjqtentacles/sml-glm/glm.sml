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
    fun mul ((a00,a10,a20, a01,a11,a21, a02,a12,a22) : t,
             (b00,b10,b20, b01,b11,b21, b02,b12,b22) : t) : t =
      let
        (* result row r, col c = sum_k A[r,k] B[k,c] *)
        val r00 = a00*b00 + a01*b10 + a02*b20
        val r10 = a10*b00 + a11*b10 + a12*b20
        val r20 = a20*b00 + a21*b10 + a22*b20
        val r01 = a00*b01 + a01*b11 + a02*b21
        val r11 = a10*b01 + a11*b11 + a12*b21
        val r21 = a20*b01 + a21*b11 + a22*b21
        val r02 = a00*b02 + a01*b12 + a02*b22
        val r12 = a10*b02 + a11*b12 + a12*b22
        val r22 = a20*b02 + a21*b12 + a22*b22
      in
        (r00,r10,r20, r01,r11,r21, r02,r12,r22)
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
            (* cofactor / adjugate, then transpose -> stored column-major *)
            val c00 =  (m11*m22 - m12*m21)
            val c01 = ~(m10*m22 - m12*m20)
            val c02 =  (m10*m21 - m11*m20)
            val c10 = ~(m01*m22 - m02*m21)
            val c11 =  (m00*m22 - m02*m20)
            val c12 = ~(m00*m21 - m01*m20)
            val c20 =  (m01*m12 - m02*m11)
            val c21 = ~(m00*m12 - m02*m10)
            val c22 =  (m00*m11 - m01*m10)
            (* inverse[r,c] = inv * cofactor[c,r] *)
            val i00 = inv*c00  val i01 = inv*c10  val i02 = inv*c20
            val i10 = inv*c01  val i11 = inv*c11  val i12 = inv*c21
            val i20 = inv*c02  val i21 = inv*c12  val i22 = inv*c22
          in
            SOME (i00,i10,i20, i01,i11,i21, i02,i12,i22)
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
            (* inverse stored column-major: index c*4+r *)
            val b00 = ( e(1,1)*c5 - e(1,2)*c4 + e(1,3)*c3) * invd
            val b01 = (ne(0,1)*c5 + e(0,2)*c4 - e(0,3)*c3) * invd
            val b02 = ( e(3,1)*s5 - e(3,2)*s4 + e(3,3)*s3) * invd
            val b03 = (ne(2,1)*s5 + e(2,2)*s4 - e(2,3)*s3) * invd
            val b10 = (ne(1,0)*c5 + e(1,2)*c2 - e(1,3)*c1) * invd
            val b11 = ( e(0,0)*c5 - e(0,2)*c2 + e(0,3)*c1) * invd
            val b12 = (ne(3,0)*s5 + e(3,2)*s2 - e(3,3)*s1) * invd
            val b13 = ( e(2,0)*s5 - e(2,2)*s2 + e(2,3)*s1) * invd
            val b20 = ( e(1,0)*c4 - e(1,1)*c2 + e(1,3)*c0) * invd
            val b21 = (ne(0,0)*c4 + e(0,1)*c2 - e(0,3)*c0) * invd
            val b22 = ( e(3,0)*s4 - e(3,1)*s2 + e(3,3)*s0) * invd
            val b23 = (ne(2,0)*s4 + e(2,1)*s2 - e(2,3)*s0) * invd
            val b30 = (ne(1,0)*c3 + e(1,1)*c1 - e(1,2)*c0) * invd
            val b31 = ( e(0,0)*c3 - e(0,1)*c1 + e(0,2)*c0) * invd
            val b32 = (ne(3,0)*s3 + e(3,1)*s1 - e(3,2)*s0) * invd
            val b33 = ( e(2,0)*s3 - e(2,1)*s1 + e(2,2)*s0) * invd
          in
            (* fromRows because b** are addressed [row,col] *)
            SOME (fromRows ((b00,b01,b02,b03),(b10,b11,b12,b13),
                            (b20,b21,b22,b23),(b30,b31,b32,b33)))
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

    fun perspective {fovy, aspect, near, far} =
      let
        val f = 1.0 / Math.tan (fovy / 2.0)
        val nf = 1.0 / (near - far)
      in
        fromRows
          (( f/aspect, 0.0, 0.0,              0.0),
           ( 0.0,      f,   0.0,              0.0),
           ( 0.0,      0.0, (far+near)*nf,    2.0*far*near*nf),
           ( 0.0,      0.0, ~1.0,             0.0))
      end

    fun ortho {left, right, bottom, top, near, far} =
      let
        val rl = 1.0 / (right - left)
        val tb = 1.0 / (top - bottom)
        val fn_ = 1.0 / (far - near)
      in
        fromRows
          (( 2.0*rl, 0.0,    0.0,      ~(right+left)*rl),
           ( 0.0,    2.0*tb, 0.0,      ~(top+bottom)*tb),
           ( 0.0,    0.0,    ~2.0*fn_, ~(far+near)*fn_),
           ( 0.0,    0.0,    0.0,      1.0))
      end

    fun lookAt {eye, center, up} =
      let
        val f = Vec3.normalize (Vec3.sub (center, eye))   (* forward *)
        val s = Vec3.normalize (Vec3.cross (f, up))       (* right *)
        val u = Vec3.cross (s, f)                         (* true up *)
        val (sx,sy,sz) = s
        val (ux,uy,uz) = u
        val (fx,fy,fz) = f
      in
        fromRows
          (( sx,  sy,  sz,  ~(Vec3.dot (s, eye))),
           ( ux,  uy,  uz,  ~(Vec3.dot (u, eye))),
           (~fx, ~fy, ~fz,    Vec3.dot (f, eye)),
           ( 0.0, 0.0, 0.0, 1.0))
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

    fun slerp (q1, q2, t) =
      let
        val d = dot (q1, q2)
        (* take shorter arc *)
        val (q2, d) = if d < 0.0 then (scale (~1.0, q2), ~d) else (q2, d)
      in
        if d > 0.9995 then
          (* nearly parallel: linear interpolate + normalize *)
          normalize (add (q1, scale (t, sub (q2, q1))))
        else
          let
            val theta0 = Math.acos (clamp (d, ~1.0, 1.0))
            val theta = theta0 * t
            val sinTheta0 = Math.sin theta0
            val s2 = Math.sin theta / sinTheta0
            val s1 = Math.cos theta - d * s2
          in
            add (scale (s1, q1), scale (s2, q2))
          end
      end
    and sub (a, b) = add (a, scale (~1.0, b))

    fun toMat4 q =
      let
        val (qw,qx,qy,qz) = normalize q
        val xx = qx*qx val yy = qy*qy val zz = qz*qz
        val xy = qx*qy val xz = qx*qz val yz = qy*qz
        val wx = qw*qx val wy = qw*qy val wz = qw*qz
      in
        Mat4.fromRows
          (( 1.0 - 2.0*(yy+zz), 2.0*(xy - wz),     2.0*(xz + wy),     0.0),
           ( 2.0*(xy + wz),     1.0 - 2.0*(xx+zz), 2.0*(yz - wx),     0.0),
           ( 2.0*(xz - wy),     2.0*(yz + wx),     1.0 - 2.0*(xx+yy), 0.0),
           ( 0.0,               0.0,               0.0,               1.0))
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
