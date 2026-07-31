//! 平面几何的最小工具集。
//!
//! 刻意只放"布局真正用得到"的东西,不铺一套通用向量库 —— 多出来的每个函数都
//! 是一处没有判据守着的表面。

/// 平面上的一个点。
///
/// # 单位是什么
///
/// **不是埃。** 2D 结构图不是分子的比例模型:芳环画成正六边形、所有键画成等长,
/// 这两件事在真实几何里都不成立。这里的长度单位就是"一个键长",取
/// [`BOND_LEN`] = 1.0,下游按自己的画布尺度缩放即可。
///
/// 把它当埃来用会诱使人拿 2D 坐标去算距离、判断接触 —— 那个数没有物理意义。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point2 {
    /// 横坐标
    pub x: f64,
    /// 纵坐标
    pub y: f64,
}

/// 标准键长。布局里所有键都画成这个长度。
pub const BOND_LEN: f64 = 1.0;

impl Point2 {
    /// 构造。
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// 原点。
    pub const ORIGIN: Self = Self { x: 0.0, y: 0.0 };

    /// 到另一点的距离。
    pub fn dist(self, other: Self) -> f64 {
        (self - other).norm()
    }

    /// 到原点的距离。
    pub fn norm(self) -> f64 {
        self.x.hypot(self.y)
    }

    /// 单位化。零向量原样返回 —— 除以零会得到 NaN,而 NaN 坐标会一路传染到
    /// 整张图,且不会在任何一步报错。
    pub fn normalized(self) -> Self {
        let n = self.norm();
        if n < f64::EPSILON {
            self
        } else {
            Self::new(self.x / n, self.y / n)
        }
    }

    /// 绕原点旋转 `radians`。
    pub fn rotated(self, radians: f64) -> Self {
        let (s, c) = radians.sin_cos();
        Self::new(self.x * c - self.y * s, self.x * s + self.y * c)
    }

    /// 绕 `pivot` 旋转 `radians`。
    pub fn rotated_about(self, pivot: Self, radians: f64) -> Self {
        (self - pivot).rotated(radians) + pivot
    }

    /// 与 `x` 轴正方向的夹角,范围 (-π, π]。
    pub fn angle(self) -> f64 {
        self.y.atan2(self.x)
    }

    /// 叉积的 z 分量。**符号就是"在直线的哪一侧"**,顺反判定全靠它。
    pub fn cross(self, other: Self) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// 点积。
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// 沿 `through` 方向、过 `pivot` 的直线做镜像。
    ///
    /// 消冲突时"翻转一根键"就是把它一侧的子树对这根键的轴做镜像。
    pub fn mirrored(self, pivot: Self, through: Self) -> Self {
        let d = through.normalized();
        let v = self - pivot;
        // v 在轴上的投影乘 2 再减 v,就是镜像
        let proj = d * (2.0 * v.dot(d));
        proj - v + pivot
    }
}

impl std::ops::Add for Point2 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y)
    }
}

impl std::ops::Sub for Point2 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y)
    }
}

impl std::ops::Mul<f64> for Point2 {
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Self::new(self.x * k, self.y * k)
    }
}

/// `p` 在有向直线 `a→b` 的哪一侧:正=左,负=右,0=线上。
///
/// 双键顺反的几何反读就是问两个参照原子是否同号。
pub fn side_of(a: Point2, b: Point2, p: Point2) -> f64 {
    (b - a).cross(p - a)
}

/// 两条线段是否**真交叉**(在内部相交)。
///
/// 共端点不算交叉 —— 相邻的两根键必然共一个端点,那不是缺陷。端点落在另一条
/// 线段内部(退化情形)也算交叉:那在图上同样是一处读不清的地方。
pub fn segments_cross(p1: Point2, p2: Point2, q1: Point2, q2: Point2) -> bool {
    // 共端点直接放过
    const EPS: f64 = 1e-9;
    for a in [p1, p2] {
        for b in [q1, q2] {
            if a.dist(b) < EPS {
                return false;
            }
        }
    }
    let d1 = side_of(q1, q2, p1);
    let d2 = side_of(q1, q2, p2);
    let d3 = side_of(p1, p2, q1);
    let d4 = side_of(p1, p2, q2);
    // 严格异号 = 真交叉
    if ((d1 > EPS && d2 < -EPS) || (d1 < -EPS && d2 > EPS))
        && ((d3 > EPS && d4 < -EPS) || (d3 < -EPS && d4 > EPS))
    {
        return true;
    }
    // 退化:某个端点落在另一条线段内部
    let on = |a: Point2, b: Point2, p: Point2, d: f64| {
        d.abs() < EPS
            && p.x >= a.x.min(b.x) - EPS
            && p.x <= a.x.max(b.x) + EPS
            && p.y >= a.y.min(b.y) - EPS
            && p.y <= a.y.max(b.y) + EPS
    };
    on(q1, q2, p1, d1) || on(q1, q2, p2, d2) || on(p1, p2, q1, d3) || on(p1, p2, q2, d4)
}

/// `p` 落在多边形 `poly` 里面吗。顶点按环上顺序给。
///
/// # 为什么不能拿"和形心同一侧"顶替
///
/// 凸多边形上两者一致,凹的就不一致了 —— 而**环并不总是凸的**:退化布局
/// (桥环)给出的环、卟啉那种大环,形心甚至可能落在环外。那时"偏向形心"会
/// 把双键的第二条线推到环外去,正是要防的那种错。实测语料里有 45 处两者判反。
///
/// 用的是奇偶规则(射线法):从 `p` 向 +x 打一条射线,数它穿过多少条边。
/// 自交的多边形也给得出答案 —— 退化布局确实会给出自交的环。
///
/// `p` 正好落在边上时结果没有定义。调用方应当从边上挪开一点再问。
pub fn point_in_polygon(p: Point2, poly: &[Point2]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        // 只数跨过 p 那条水平线的边。用 `>` 的异或写法把顶点正好在线上的
        // 情形算作"属于上面那条边",于是每个顶点只被数一次。
        if (a.y > p.y) != (b.y > p.y) {
            let t = (p.y - a.y) / (b.y - a.y);
            if a.x + t * (b.x - a.x) > p.x {
                inside = !inside;
            }
        }
    }
    inside
}

/// `n` 边正多边形的顶点,边长为 [`BOND_LEN`],中心在原点。
///
/// `start_angle` 是第 0 个顶点相对 x 轴的角度。
pub fn regular_polygon(n: usize, start_angle: f64) -> Vec<Point2> {
    debug_assert!(n >= 3, "多边形至少三个顶点");
    // 边长 s 与外接圆半径 r 的关系:s = 2 r sin(π/n)
    let r = BOND_LEN / (2.0 * (std::f64::consts::PI / n as f64).sin());
    let step = std::f64::consts::TAU / n as f64;
    (0..n)
        .map(|i| Point2::new(r, 0.0).rotated(start_angle + step * i as f64))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    #[test]
    fn a_regular_polygon_has_unit_sides() {
        // 正多边形的边长必须**正好是一个键长** —— 这是"键长全等"那条判据在
        // 环上的来源。半径公式写错(比如把 sin 写成 tan)不会让程序出错,
        // 只会让每个环都稍微大一点或小一点,肉眼很难发现。
        for n in 3..=12 {
            let p = regular_polygon(n, 0.0);
            assert_eq!(p.len(), n);
            for i in 0..n {
                let d = p[i].dist(p[(i + 1) % n]);
                assert!((d - BOND_LEN).abs() < TOL, "{n} 边形第 {i} 条边长 {d}");
            }
        }
    }

    #[test]
    fn mirroring_twice_is_identity() {
        let pivot = Point2::new(1.0, 2.0);
        let axis = Point2::new(0.6, -0.8);
        let p = Point2::new(-3.0, 5.0);
        let back = p.mirrored(pivot, axis).mirrored(pivot, axis);
        assert!(p.dist(back) < TOL);
    }

    #[test]
    fn mirroring_flips_which_side_a_point_is_on() {
        // 镜像必须真的把点换到另一侧。写成"投影"而不是"两倍投影减自身"
        // 会得到一个把点压到轴上的操作 —— 它同样不报错,但翻转失效,
        // 消冲突就会静默地什么也不做。
        let a = Point2::ORIGIN;
        let b = Point2::new(1.0, 0.0);
        let p = Point2::new(0.5, 1.0);
        let m = p.mirrored(a, b - a);
        assert!(
            side_of(a, b, p) * side_of(a, b, m) < 0.0,
            "镜像后没有换侧:{m:?}"
        );
    }

    #[test]
    fn crossing_segments_are_detected_and_touching_ones_are_not() {
        let cross = segments_cross(
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 2.0),
            Point2::new(0.0, 2.0),
            Point2::new(2.0, 0.0),
        );
        assert!(cross, "对角线必须判为交叉");

        // 共端点的两根键(化学里相邻的键)不算交叉
        let shared = segments_cross(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.5, 1.0),
        );
        assert!(!shared, "共端点不该判为交叉");

        let apart = segments_cross(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(1.0, 1.0),
        );
        assert!(!apart, "平行且分离的两段不该判为交叉");
    }

    #[test]
    fn a_point_inside_a_convex_polygon_is_found_and_one_outside_is_not() {
        let hex = regular_polygon(6, 0.0);
        assert!(point_in_polygon(Point2::ORIGIN, &hex), "形心该在里面");
        assert!(
            !point_in_polygon(Point2::new(5.0, 0.0), &hex),
            "远处的点该在外面"
        );
        // 同一条水平线上、多边形左边的点:射线法数右侧交点,这里是 2 个 → 外面
        assert!(
            !point_in_polygon(Point2::new(-5.0, 0.0), &hex),
            "左边也是外面"
        );
    }

    #[test]
    fn the_centroid_of_a_concave_polygon_can_be_outside_it() {
        // **这条才是要这个函数的理由。** 凸多边形上"在里面"和"和形心同侧"
        // 一致,凹的就不一致 —— 而环并不总是凸的。
        //
        // 一个 U 形:形心落在 U 的口子里,也就是多边形**外面**。
        let u = [
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 3.0),
            Point2::new(2.0, 3.0),
            Point2::new(2.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 3.0),
            Point2::new(0.0, 3.0),
        ];
        let c = u.iter().fold(Point2::ORIGIN, |s, p| s + *p) * (1.0 / u.len() as f64);
        assert!(
            !point_in_polygon(c, &u),
            "这个 U 形的形心 {c:?} 落在口子里,该判成外面"
        );
        assert!(point_in_polygon(Point2::new(0.5, 2.0), &u), "左臂里面");
        assert!(point_in_polygon(Point2::new(1.5, 0.5), &u), "底下里面");
        assert!(!point_in_polygon(Point2::new(1.5, 2.0), &u), "口子里是外面");
    }

    #[test]
    fn a_degenerate_polygon_is_refused_instead_of_panicking() {
        // 少于三个顶点构不成多边形。退化布局里什么都可能传进来,不能 panic。
        assert!(!point_in_polygon(Point2::ORIGIN, &[]));
        assert!(!point_in_polygon(
            Point2::ORIGIN,
            &[Point2::new(1.0, 1.0), Point2::new(2.0, 2.0)]
        ));
    }

    #[test]
    fn normalizing_a_zero_vector_does_not_produce_nan() {
        // NaN 坐标会一路传染到整张图,而且沿途一步都不会报错
        let z = Point2::ORIGIN.normalized();
        assert!(z.x.is_finite() && z.y.is_finite());
    }
}
