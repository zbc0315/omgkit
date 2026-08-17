//! 三维向量与"按内坐标摆一个原子"。
//!
//! 这个模块只做几何,不认识分子。

/// 三维点 / 向量。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    /// x
    pub x: f64,
    /// y
    pub y: f64,
    /// z
    pub z: f64,
}

impl Point3 {
    /// 原点。
    pub const ORIGIN: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// 造一个点。
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// 长度。
    #[must_use]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// 点积。
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x.mul_add(o.x, self.y.mul_add(o.y, self.z * o.z))
    }

    /// 叉积。
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y.mul_add(o.z, -(self.z * o.y)),
            self.z.mul_add(o.x, -(self.x * o.z)),
            self.x.mul_add(o.y, -(self.y * o.x)),
        )
    }

    /// 归一化。**长度为 0 时返回 `None`**,不返回 NaN。
    #[must_use]
    pub fn normalized(self) -> Option<Self> {
        let n = self.norm();
        if n < 1e-12 || !n.is_finite() {
            return None;
        }
        Some(self * (1.0 / n))
    }

    /// 两点距离。
    #[must_use]
    pub fn dist(self, o: Self) -> f64 {
        (self - o).norm()
    }

    /// 三个坐标都是有限数。
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
}

impl std::ops::Add for Point3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl std::ops::Sub for Point3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl std::ops::Mul<f64> for Point3 {
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Self::new(self.x * k, self.y * k, self.z * k)
    }
}

/// `a`–`b`–`c` 的夹角(弧度,顶点在 `b`)。三点退化时返回 `None`。
#[must_use]
pub fn angle_at(a: Point3, b: Point3, c: Point3) -> Option<f64> {
    let u = (a - b).normalized()?;
    let v = (c - b).normalized()?;
    Some(u.dot(v).clamp(-1.0, 1.0).acos())
}

/// `a`–`b`–`c`–`d` 的二面角(弧度,范围 `(−π, π]`)。退化时返回 `None`。
#[must_use]
pub fn dihedral(a: Point3, b: Point3, c: Point3, d: Point3) -> Option<f64> {
    // **符号是 IUPAC 的那一支。** 头一版把 `b1` 取成 `b − a`、用 `n1 × ê` 当 y 基,
    // 结果整体**差一个负号**(手算 +90° 的例子它给 −90°)—— 而 NeRF 那边是对的,
    // 差点被我改错。所以下面这条与 `place_nerf` 是配对的,**改一个就要改另一个**,
    // 判据是 `dihedral_matches_a_hand_computed_case`。
    let e = (c - b).normalized()?;
    let b0 = a - b;
    let b3 = d - c;
    let v = b0 - e * b0.dot(e);
    let w = b3 - e * b3.dot(e);
    let x = v.dot(w);
    let y = e.cross(v).dot(w);
    if x.abs() < 1e-30 && y.abs() < 1e-30 {
        return None;
    }
    Some(y.atan2(x))
}

/// **按内坐标摆第四个原子**(NeRF,Natural Extension Reference Frame)。
///
/// 已知 `a`、`b`、`c` 三个点,要摆的新点 `d` 满足:
/// `|cd| = bond`、`∠bcd = angle`、二面角 `a-b-c-d = torsion`(都是弧度)。
///
/// # 为什么用这个而不是解方程
///
/// 它是**闭式**的:一次三角函数 + 一次标架变换,没有迭代、没有收敛判据、
/// 不会失败(除非输入本身退化)。链上的每个原子都靠它摆,`O(N)`。
///
/// # 退化
///
/// `a`、`b`、`c` 共线时定不出标架,返回 `None` —— 调用方要自己兜底
/// (通常是换一个参考原子,或对第一个原子用固定方向)。
#[must_use]
pub fn place_nerf(
    a: Point3,
    b: Point3,
    c: Point3,
    bond: f64,
    angle: f64,
    torsion: f64,
) -> Option<Point3> {
    if !(bond.is_finite() && angle.is_finite() && torsion.is_finite()) || bond <= 0.0 {
        return None;
    }
    // 局部坐标系:x 沿 c←b,z 垂直于 abc 平面
    let bc = (c - b).normalized()?;
    let ab = b - a;
    let n = ab.cross(bc).normalized()?; // 共线时这里 None
                                        // **中间那个基是 `n × bc`,不是 `bc × n`** —— 写反了二面角整体差 180°
                                        // (单元测试第一次跑就逮到:要 −179° 给出 +1°)。
    let m = n.cross(bc);
    // 新点在局部系里的方向
    let (sa, ca) = angle.sin_cos();
    let (st, ct) = torsion.sin_cos();
    let d2 = Point3::new(-ca, sa * ct, sa * st);
    let dir = bc * d2.x + m * d2.y + n * d2.z;
    let p = c + dir * bond;
    p.is_finite().then_some(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// NeRF 摆出来的点必须**精确**满足给的键长、键角、二面角。
    ///
    /// 这是整条链的地基:它错了后面全错,而且错得很隐蔽(几何看着还像分子)。
    #[test]
    fn nerf_reproduces_the_internal_coordinates_it_was_given() {
        let a = Point3::new(1.0, 0.4, -0.2);
        let b = Point3::new(0.0, 0.0, 0.0);
        let c = Point3::new(1.5, 0.0, 0.0);
        let mut checked = 0;
        for &bond in &[1.09, 1.4, 1.54] {
            for &deg in &[60.0, 90.0, 109.471, 120.0, 150.0] {
                for &tor in &[-179.0, -90.0, 0.0, 60.0, 120.0, 180.0] {
                    let (ang, t) = (f64::to_radians(deg), f64::to_radians(tor));
                    let d = place_nerf(a, b, c, bond, ang, t).expect("该摆得出来");
                    assert!(
                        (c.dist(d) - bond).abs() < 1e-12,
                        "键长 {} 该是 {bond}",
                        c.dist(d)
                    );
                    let got_a = angle_at(b, c, d).expect("角该算得出");
                    assert!(
                        (got_a - ang).abs() < 1e-12,
                        "键角 {:.6}° 该是 {deg}°",
                        got_a.to_degrees()
                    );
                    let got_t = dihedral(a, b, c, d).expect("二面角该算得出");
                    // 二面角要按**绕回**比:|Δ| 与 2π−|Δ| 取小
                    let raw = (got_t - t).abs();
                    let diff = raw.min(std::f64::consts::TAU - raw);
                    assert!(diff < 1e-9, "二面角 {:.6}° 该是 {tor}°", got_t.to_degrees());
                    checked += 1;
                }
            }
        }
        assert!(checked >= 80, "只验了 {checked} 组 —— 这条测试快变恒真了");
    }

    /// **先把量角的尺子钉死,再拿它去验 NeRF。**
    ///
    /// 手算的例子:`a=(0,1,0)`、`b=原点`、`c=(1,0,0)`、`d=(1,0,1)`。
    /// 沿 `b→c`(+x)看过去,`a` 的投影在 +y、`d` 的投影在 +z,
    /// +y → +z 绕 +x 是**右手 +90°**。
    ///
    /// 这条是必须的:头一版 `dihedral` 的符号是反的(它给 −90°),
    /// 而我差点因此去改本来正确的 `place_nerf`。
    #[test]
    fn dihedral_matches_a_hand_computed_case() {
        let a = Point3::new(0.0, 1.0, 0.0);
        let b = Point3::ORIGIN;
        let c = Point3::new(1.0, 0.0, 0.0);
        let d = Point3::new(1.0, 0.0, 1.0);
        let got = dihedral(a, b, c, d).expect("该算得出").to_degrees();
        assert!((got - 90.0).abs() < 1e-9, "该是 +90°,得到 {got:.6}°");
        // 把 d 镜像到 −z,应当变成 −90°
        let d2 = Point3::new(1.0, 0.0, -1.0);
        let got2 = dihedral(a, b, c, d2).expect("该算得出").to_degrees();
        assert!((got2 + 90.0).abs() < 1e-9, "该是 −90°,得到 {got2:.6}°");
    }

    /// 三点共线时定不出标架,必须说定不出,**不能返回 NaN**。
    #[test]
    fn nerf_says_no_when_the_reference_atoms_are_collinear() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(2.0, 0.0, 0.0);
        assert!(place_nerf(a, b, c, 1.5, 2.0, 1.0).is_none());
    }

    /// 反式(180°)的丁烷骨架:四个碳该在一个平面上,且首尾最远。
    #[test]
    fn an_anti_backbone_is_planar_and_extended() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.54, 0.0, 0.0);
        let ang = f64::to_radians(109.471);
        let c = place_nerf(Point3::new(0.0, 1.0, 0.0), a, b, 1.54, ang, 0.0).expect("c");
        let d = place_nerf(a, b, c, 1.54, ang, std::f64::consts::PI).expect("d");
        let t = dihedral(a, b, c, d).expect("二面角");
        assert!(
            (t.abs() - std::f64::consts::PI).abs() < 1e-9,
            "该是 180°,得到 {:.3}°",
            t.to_degrees()
        );
        // **这个数是算出来的,不是猜的。** 全反式锯齿:沿链轴每段 `L·sin(θ/2)`,
        // 垂直半幅 `L·cos(θ/2)/2`(**是一半** —— 头一版把它当成整段位移,
        // 算出 4.170 而 NeRF 给 3.876,差点去改本来正确的 NeRF)。
        // 于是 C1–C4 = √((3L·sin(θ/2))² + (L·cos(θ/2))²) = 3.8756。
        let anti = a.dist(d);
        assert!(
            (anti - 3.8756).abs() < 1e-3,
            "反式首尾该是 3.8756,得到 {anti:.4}"
        );
        // 与邻位交叉(60°,2.9489)拉开 —— 不然这条判据分不出扭转角对不对
        let g = place_nerf(a, b, c, 1.54, ang, f64::to_radians(60.0)).expect("gauche");
        assert!(
            a.dist(g) < anti - 0.8,
            "邻位交叉该明显更短:{:.4} vs 反式 {anti:.4}",
            a.dist(g)
        );
    }
}
