(* glm.sig

   Pure linear algebra and transforms for graphics, in the spirit of GLM/GLSL.

   Conventions:
   - Matrices are stored column-major, matching OpenGL's memory layout, and
     transforms follow the right-handed coordinate system. Angles are in
     radians; use `radians`/`degrees` to convert.
   - `Mat*.mulV` applies a matrix to a (column) vector: result = M * v.
   - `inverse` returns NONE for singular matrices.
   - `normalize` of a zero-length vector returns the zero vector (it never
     divides by zero).

   All operations are pure and total; arithmetic is over `real`. *)

signature GLM =
sig
  val pi : real
  val radians : real -> real            (* degrees -> radians *)
  val degrees : real -> real            (* radians -> degrees *)
  val clamp   : real * real * real -> real   (* clamp (x, lo, hi) *)
  val lerp    : real * real * real -> real   (* lerp (a, b, t) *)

  structure Vec2 :
  sig
    type t
    val v : real * real -> t
    val x : t -> real
    val y : t -> real
    val zero : t
    val toList : t -> real list

    val add   : t * t -> t
    val sub   : t * t -> t
    val mulc  : t * t -> t            (* componentwise (Hadamard) product *)
    val scale : real * t -> t
    val neg   : t -> t

    val dot      : t * t -> real
    val length   : t -> real
    val lengthSq : t -> real
    val dist     : t * t -> real
    val normalize : t -> t
    val lerp : t * t * real -> t

    val equal  : t * t -> bool
    val approx : real -> t * t -> bool
    val toString : t -> string
  end

  structure Vec3 :
  sig
    type t
    val v : real * real * real -> t
    val x : t -> real
    val y : t -> real
    val z : t -> real
    val zero : t
    val toList : t -> real list

    val add   : t * t -> t
    val sub   : t * t -> t
    val mulc  : t * t -> t
    val scale : real * t -> t
    val neg   : t -> t

    val dot      : t * t -> real
    val cross    : t * t -> t
    val length   : t -> real
    val lengthSq : t -> real
    val dist     : t * t -> real
    val normalize : t -> t
    val lerp : t * t * real -> t

    val equal  : t * t -> bool
    val approx : real -> t * t -> bool
    val toString : t -> string
  end

  structure Vec4 :
  sig
    type t
    val v : real * real * real * real -> t
    val x : t -> real
    val y : t -> real
    val z : t -> real
    val w : t -> real
    val zero : t
    val toList : t -> real list

    val add   : t * t -> t
    val sub   : t * t -> t
    val mulc  : t * t -> t
    val scale : real * t -> t
    val neg   : t -> t

    val dot      : t * t -> real
    val length   : t -> real
    val lengthSq : t -> real
    val dist     : t * t -> real
    val normalize : t -> t
    val lerp : t * t * real -> t

    val equal  : t * t -> bool
    val approx : real -> t * t -> bool
    val toString : t -> string
  end

  structure Mat2 :
  sig
    type t
    val id : t
    (* Build from rows (the natural reading order); stored column-major. *)
    val fromRows : Vec2.t * Vec2.t -> t
    val fromCols : Vec2.t * Vec2.t -> t
    (* Column-major 4-element list. *)
    val fromList : real list -> t
    val toList   : t -> real list

    val add : t * t -> t
    val sub : t * t -> t
    val mul : t * t -> t
    val scale : real * t -> t
    val transpose : t -> t
    val det : t -> real
    val inverse : t -> t option
    val mulV : t * Vec2.t -> Vec2.t

    val equal  : t * t -> bool
    val approx : real -> t * t -> bool
    val toString : t -> string
  end

  structure Mat3 :
  sig
    type t
    val id : t
    (* Build from rows (the natural reading order); stored column-major. *)
    val fromRows : Vec3.t * Vec3.t * Vec3.t -> t
    val fromCols : Vec3.t * Vec3.t * Vec3.t -> t
    (* Column-major 9-element list. *)
    val fromList : real list -> t
    val toList   : t -> real list

    val add : t * t -> t
    val sub : t * t -> t
    val mul : t * t -> t
    val scale : real * t -> t
    val transpose : t -> t
    val det : t -> real
    val inverse : t -> t option
    val mulV : t * Vec3.t -> Vec3.t

    val equal  : t * t -> bool
    val approx : real -> t * t -> bool
    val toString : t -> string
  end

  structure Mat4 :
  sig
    type t
    val id : t
    val fromRows : Vec4.t * Vec4.t * Vec4.t * Vec4.t -> t
    val fromCols : Vec4.t * Vec4.t * Vec4.t * Vec4.t -> t
    val fromList : real list -> t      (* column-major, 16 elements *)
    val toList   : t -> real list

    val add : t * t -> t
    val sub : t * t -> t
    val mul : t * t -> t
    val scale : real * t -> t
    val transpose : t -> t
    val det : t -> real
    val inverse : t -> t option
    val mulV : t * Vec4.t -> Vec4.t
    (* Apply to a point (w=1) / direction (w=0), returning a Vec3. *)
    val transformPoint : t * Vec3.t -> Vec3.t
    val transformDir   : t * Vec3.t -> Vec3.t

    (* Affine builders. *)
    val translate : Vec3.t -> t
    val scaleM    : Vec3.t -> t
    val rotate    : real * Vec3.t -> t     (* angle (radians), axis *)
    val rotateX   : real -> t
    val rotateY   : real -> t
    val rotateZ   : real -> t

    (* Projection / view (right-handed, GL clip space z in [-1, 1]). *)
    val perspective : { fovy : real, aspect : real, near : real, far : real } -> t
    val ortho : { left : real, right : real, bottom : real, top : real,
                  near : real, far : real } -> t
    (* General perspective frustum (`glFrustum`): right-handed, GL clip-space
       z in [-1, 1]. `perspective` is the symmetric special case. *)
    val frustum : { left : real, right : real, bottom : real, top : real,
                    near : real, far : real } -> t
    val lookAt : { eye : Vec3.t, center : Vec3.t, up : Vec3.t } -> t

    val equal  : t * t -> bool
    val approx : real -> t * t -> bool
    val toString : t -> string
  end

  structure Quat :
  sig
    type t
    (* Components in (w, x, y, z) order. *)
    val quat : real * real * real * real -> t
    val id : t
    val w : t -> real
    val x : t -> real
    val y : t -> real
    val z : t -> real
    val toList : t -> real list

    val fromAxisAngle : Vec3.t * real -> t   (* axis, angle (radians) *)
    val add   : t * t -> t
    val scale : real * t -> t
    val mul   : t * t -> t
    val conj  : t -> t
    val dot   : t * t -> real
    val length : t -> real
    val normalize : t -> t
    val rotateV : t * Vec3.t -> Vec3.t
    val slerp : t * t * real -> t
    val toMat4  : t -> Mat4.t
    val fromMat3 : Mat3.t -> t

    val equal  : t * t -> bool
    val approx : real -> t * t -> bool
    val toString : t -> string
  end

  (* Map an object-space point to window coordinates given the model(view) and
     projection matrices and a viewport (`gluProject`): window x/y are in
     pixels, window z (depth) is in [0, 1]. *)
  val project :
    { obj : Vec3.t, model : Mat4.t, proj : Mat4.t,
      viewport : { x : real, y : real, width : real, height : real } }
    -> Vec3.t

  (* Inverse of `project` (`gluUnProject`): map a window coordinate (with z in
     [0, 1]) back to object space. Returns `Vec3.zero` if `proj * model` is
     singular or the result is at infinity, so it is total. *)
  val unproject :
    { win : Vec3.t, model : Mat4.t, proj : Mat4.t,
      viewport : { x : real, y : real, width : real, height : real } }
    -> Vec3.t
end
