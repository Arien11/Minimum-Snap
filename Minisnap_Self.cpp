#include <iostream>
#include <traj_utils/polynomial_traj.h>

PolynomialTraj PolynomialTraj::minSnapTraj(const Eigen::MatrixXd &Pos, const Eigen::Vector3d &start_vel,
                                           const Eigen::Vector3d &end_vel, const Eigen::Vector3d &start_acc,
                                           const Eigen::Vector3d &end_acc, const Eigen::VectorXd &Time)
{
  /* 约束
    轨迹起点速度、加速度
    轨迹终点速度、加速度
    轨迹中每个端点的位置连续性P_i(T) = P_(i+1)(0)
    输入的Pos:(3(三轴), 分段点的位置信息)
  */

  int seg_num = Time.size();
  Eigen::MatrixXd poly_coeff(seg_num, 3 * 6);   // 6项，那不是只是minJerk吗
  Eigen::VectorXd Px(6 * seg_num), Py(6 * seg_num), Pz(6 * seg_num);

  int num_f, num_p; // f表示定值向量, p表示优化向量
  int num_d;        // number of all segments' derivatives d表示

  // 内嵌的阶乘函数定义
  const static auto Factorial = [](int x) {
    int fac = 1;
    for (int i = x; i > 0; i--)
      fac = fac * i;
    return fac;
  };

  /* ---------- end point derivative ---------- */
  // 存放约束的向量 Dx: seg_num * (当前点的位置, 下个点的位置, 0, 0, 0, 0) 除了起点和终点有速度加速度约束其他点均为0
  Eigen::VectorXd Dx = Eigen::VectorXd::Zero(seg_num * 6);    
  Eigen::VectorXd Dy = Eigen::VectorXd::Zero(seg_num * 6);
  Eigen::VectorXd Dz = Eigen::VectorXd::Zero(seg_num * 6);

  // 存入连续性约束与微分约束到Dx，Dy，Dz
  for (int k = 0; k < seg_num; k++)
  {
    /* position to derivative */
    // 端点位置约束及连续性约束(P_0(T) = P_1(0))
    // Dx(0) = Pos(0, 0), Dx(1) = Pos(0, 1)
    // Dx(6) = Pos(0, 1), Dx(7) = Pos(0, 2)
    Dx(k * 6) = Pos(0, k);              
    Dx(k * 6 + 1) = Pos(0, k + 1);  
    Dy(k * 6) = Pos(1, k);
    Dy(k * 6 + 1) = Pos(1, k + 1);
    Dz(k * 6) = Pos(2, k);
    Dz(k * 6 + 1) = Pos(2, k + 1);

    // 起/终点的速度与加速度约束
    if (k == 0)
    {
      Dx(k * 6 + 2) = start_vel(0); // 起点的速度约束     
      Dy(k * 6 + 2) = start_vel(1);
      Dz(k * 6 + 2) = start_vel(2);

      Dx(k * 6 + 4) = start_acc(0); // 起点的加速度约束 
      Dy(k * 6 + 4) = start_acc(1);
      Dz(k * 6 + 4) = start_acc(2);
    }
    else if (k == seg_num - 1)
    {
      Dx(k * 6 + 3) = end_vel(0);  // 终点的加速度约束 
      Dy(k * 6 + 3) = end_vel(1);
      Dz(k * 6 + 3) = end_vel(2);

      Dx(k * 6 + 5) = end_acc(0);  // 终点的加速度约束 
      Dy(k * 6 + 5) = end_acc(1);
      Dz(k * 6 + 5) = end_acc(2);
    }
  }

  /* ---------- Mapping Matrix A 处理微分约束 ---------- */
  // 多项式系数Poly_Coeff映射到约束量
  Eigen::MatrixXd Ab;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);
  // 工程上是按照阶次进行排布，而非点约束顺序排布
  // 因此0,1行对应的是0阶导数，2,3行为1阶导数，4,5行为2阶导数(就是起点和终点的多项式导数作为两行)
  // 由于微分约束针对是起点与终点，起点为t=0(t=0代到多项式中很多都没了)，则对应的每行有效值只有导完之后t阶次为0的那项，呈阶梯状排布
  for (int k = 0; k < seg_num; k++)   // 这个就是A
  {
    Ab = Eigen::MatrixXd::Zero(6, 6);   // 这个就是Am
    // 是针对jerk(位置的三阶导)的，故而只做到位置的二阶导，即最多到加速度的约束
    for (int i = 0; i < 3; i++)         // i对应当前为几阶导
    {
      // 多项式的每阶导数对应一个约束
      Ab(2 * i, i) = Factorial(i);      // 
      for (int j = i; j < 6; j++)    // 多项式为5阶多项式，根据i阶导进行计算导数系数
        Ab(2 * i + 1, j) = Factorial(j) / Factorial(j - i) * pow(Time(k), j - i);   
    }
    A.block(k * 6, k * 6, 6, 6) = Ab;   // 将第k段的映射矩阵赋值给大矩阵A

    // Ab(0, 0) = Factorial(0)， 起点位置约束
    // Ab(2, 1) = Factorial(1)， 起点速度约束
    // Ab(4, 2) = Factorial(2)   起点加速度约束  

    // Ab(1, j) = 1* t_k^(j-0), j=i,...,5    终点位置约束
    // Ab(3, j) = Factorial(j) / Factorial(j - i) * t_k^(j-1), j=i,...,5 终点速度约束
    // Ab(3, j) = Factorial(j) / Factorial(j - i) * t_k^(j-2), j=i,...,5 终点加速度约束

  }

  /* ---------- Produce Selection Matrix C 处理连续性约束 ---------- */
  Eigen::MatrixXd Ct, C;
  // 置换矩阵C
  num_f = 2 * seg_num + 4; // 3 + 3 + (seg_num - 1) * 2 = 2m + 4     每段轨迹中两个端点的位置约束，即(s-1) * seg_num个, 加上起点终点的速度、加速度约束，即4个
  num_p = 2 * seg_num - 2; // (seg_num - 1) * 2 = 2m - 2，           每段轨迹中(除开起点和终点)两个端点的待优化约束，(s-1) * (seg_num-1)
  num_d = 6 * seg_num;     // 2x3 + 6 * (seg_num - 1) = 6 * seg_num  起/终点约束，中间点左右连续性约束
  Ct = Eigen::MatrixXd::Zero(num_d, num_f + num_p);     

  // d:  各个段的按导的阶数排，阶数内按时间顺序排
  /*
      d = [d_(0,0)^(0), d(0, T)^(0), d_(0,0)^(1), d(0, T)^(1), d_(0,0)^(2), d(0, T)^(2), d_(1,0)^(0), d(1, T)^(0),d_(1,0)^(1), d(1, T)^(1)...]^T
      
      df = [d_(0,0)^(0), d_(0,0)^(1), d_(0,0)^(2), d(0, T)^(0), 

            d(1,0)^(0),  d(1,T)^(0), 
            d(2,0)^(0),  d(2,T)^(0),
                      ...
            d(seg_num-2,0)^(0), d_(seg_num-2, T)^(0),   

            d(seg_num-1, 0)^(0), d(seg_num-1, T)^(0), d(seg_num-1, T)^(1), d(seg_num-1, T)^(2)]^T

      dp = [d(1,0)^(1), d(1,0)^(2), d(1,T)^(1), d(1,T)^(2), 
            d(2,0)^(1), d(2,0)^(2), d(2,T)^(1), d(2,T)^(2),
            ...
            d(seg_num-2,0)^(1), d(seg_num-2,0)^(2), d(seg_num-2,T)^(1), d(seg_num-2,T)^(2)]^T
      
  */

  // 对端点阶数分组，阶数内按时间顺序 
  // 例子：对应0阶导的有两行，第一行对应时间t=0，第二行对应时间t=T

  // 每段有6行，每两行为一阶导
  // 每段0,2,4为t=0
  // 每段1,3,5为t=T

  // 起点微分约束，起点约束是用于t=0的,因而如下
  Ct(0, 0) = 1;
  Ct(2, 1) = 1;
  Ct(4, 2) = 1; // stack the start point

  // 起始段的末端点的连续性约束
  Ct(1, 3) = 1;                                   // 位置约束是定值
  Ct(3, 2 * seg_num + 4) = 1;                     // 速度与加速度约束在dp中
  Ct(5, 2 * seg_num + 4 + 1) = 1; 

  // 终点约束，终点约束是针对于t=T的,因而如下
  Ct(6 * (seg_num - 1) + 1, 2 * seg_num + 1) = 1; // Stack the end point
  Ct(6 * (seg_num - 1) + 3, 2 * seg_num + 2) = 1; // Stack the end point
  Ct(6 * (seg_num - 1) + 5, 2 * seg_num + 3) = 1; // Stack the end point

  Ct(6 * (seg_num - 1) + 0, 2 * seg_num + 0) = 1; // 位置约束是定值
  Ct(6 * (seg_num - 1) + 2, 4 * seg_num + 0) = 1; // 速度与加速度约束在dp中 d_(seg_num - 1, 0)^(1) = d_(seg_num - 2, T)^(1), 就dp中的倒数两个 f+p - 2 = 2 * seg_num + 4 + 2 * seg_num - 2 -2 = 4s
  Ct(6 * (seg_num - 1) + 4, 4 * seg_num + 1) = 1; // d_(seg_num - 1, 0)^(2) = d_(seg_num - 2, T)^(2)


  // 从第二段开始, j-1是因为下标从0开始的其实也就是第2段
  // df中除了开头结尾，中间的都是2个为一组的，随便找
  // dp中每段是以四个为一组(开头的速度加速度2个，结尾的速度加速度2个)
  for (int j = 2; j < seg_num; j++)
  {
    // 位置都是定值在df中，其他都是待优化的在dp中

    // t = 0
    Ct(6 * (j - 1) + 0, 2 + 2 * (j - 1) + 0) = 1;    
    Ct(6 * (j - 1) + 2, 2 * seg_num + 4 + 2 * (j - 2) + 0) = 1;
    Ct(6 * (j - 1) + 4, 2 * seg_num + 4 + 2 * (j - 2) + 1) = 1;

    // t = T
    Ct(6 * (j - 1) + 1, 2 + 2 * (j - 1) + 1) = 1;
    Ct(6 * (j - 1) + 3, 2 * seg_num + 4 + 2 * (j - 1) + 0) = 1;
    Ct(6 * (j - 1) + 5, 2 * seg_num + 4 + 2 * (j - 1) + 1) = 1;
  }

  C = Ct.transpose();

  Eigen::VectorXd Dx1 = C * Dx;
  Eigen::VectorXd Dy1 = C * Dy;
  Eigen::VectorXd Dz1 = C * Dz;

  /* ---------- minimum snap matrix ---------- */
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);

  for (int k = 0; k < seg_num; k++)
  {
    for (int i = 3; i < 6; i++)
    {
      for (int j = 3; j < 6; j++)
      {
        Q(k * 6 + i, k * 6 + j) =
            i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) / (i + j - 5) * pow(Time(k), (i + j - 5));
      }
    }
  }

  /* ---------- R matrix ---------- */
  Eigen::MatrixXd R = C * A.transpose().inverse() * Q * A.inverse() * Ct;

  Eigen::VectorXd Dxf(2 * seg_num + 4), Dyf(2 * seg_num + 4), Dzf(2 * seg_num + 4);

  Dxf = Dx1.segment(0, 2 * seg_num + 4);
  Dyf = Dy1.segment(0, 2 * seg_num + 4);
  Dzf = Dz1.segment(0, 2 * seg_num + 4);

  Eigen::MatrixXd Rff(2 * seg_num + 4, 2 * seg_num + 4);
  Eigen::MatrixXd Rfp(2 * seg_num + 4, 2 * seg_num - 2);
  Eigen::MatrixXd Rpf(2 * seg_num - 2, 2 * seg_num + 4);
  Eigen::MatrixXd Rpp(2 * seg_num - 2, 2 * seg_num - 2);

  Rff = R.block(0, 0, 2 * seg_num + 4, 2 * seg_num + 4);
  Rfp = R.block(0, 2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2);
  Rpf = R.block(2 * seg_num + 4, 0, 2 * seg_num - 2, 2 * seg_num + 4);
  Rpp = R.block(2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2, 2 * seg_num - 2);

  /* ---------- close form solution ---------- */

  Eigen::VectorXd Dxp(2 * seg_num - 2), Dyp(2 * seg_num - 2), Dzp(2 * seg_num - 2);
  Dxp = -(Rpp.inverse() * Rfp.transpose()) * Dxf;
  Dyp = -(Rpp.inverse() * Rfp.transpose()) * Dyf;
  Dzp = -(Rpp.inverse() * Rfp.transpose()) * Dzf;

  Dx1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dxp;
  Dy1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dyp;
  Dz1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dzp;

  Px = (A.inverse() * Ct) * Dx1;
  Py = (A.inverse() * Ct) * Dy1;
  Pz = (A.inverse() * Ct) * Dz1;

  // 得到最后的多项式矩阵
  for (int i = 0; i < seg_num; i++)
  {
    poly_coeff.block(i, 0, 1, 6) = Px.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 6, 1, 6) = Py.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 12, 1, 6) = Pz.segment(i * 6, 6).transpose();
  }

  /* ---------- use polynomials ---------- */
  PolynomialTraj poly_traj;
  // 多项式中解算出轨迹
  for (int i = 0; i < poly_coeff.rows(); ++i)
  {
    vector<double> cx(6), cy(6), cz(6);
    for (int j = 0; j < 6; ++j)
    {
      cx[j] = poly_coeff(i, j), cy[j] = poly_coeff(i, j + 6), cz[j] = poly_coeff(i, j + 12);
    }
    reverse(cx.begin(), cx.end());
    reverse(cy.begin(), cy.end());
    reverse(cz.begin(), cz.end());
    double ts = Time(i);
    poly_traj.addSegment(cx, cy, cz, ts);
  }

  return poly_traj;
}

PolynomialTraj PolynomialTraj::one_segment_traj_gen(const Eigen::Vector3d &start_pt, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                                    const Eigen::Vector3d &end_pt, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc,
                                                    double t)
{
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6), Crow(1, 6);
  Eigen::VectorXd Bx(6), By(6), Bz(6);

  C(0, 5) = 1;
  C(1, 4) = 1;
  C(2, 3) = 2;
  Crow << pow(t, 5), pow(t, 4), pow(t, 3), pow(t, 2), t, 1;
  C.row(3) = Crow;
  Crow << 5 * pow(t, 4), 4 * pow(t, 3), 3 * pow(t, 2), 2 * t, 1, 0;
  C.row(4) = Crow;
  Crow << 20 * pow(t, 3), 12 * pow(t, 2), 6 * t, 2, 0, 0;
  C.row(5) = Crow;

  Bx << start_pt(0), start_vel(0), start_acc(0), end_pt(0), end_vel(0), end_acc(0);
  By << start_pt(1), start_vel(1), start_acc(1), end_pt(1), end_vel(1), end_acc(1);
  Bz << start_pt(2), start_vel(2), start_acc(2), end_pt(2), end_vel(2), end_acc(2);

  Eigen::VectorXd Cofx = C.colPivHouseholderQr().solve(Bx);
  Eigen::VectorXd Cofy = C.colPivHouseholderQr().solve(By);
  Eigen::VectorXd Cofz = C.colPivHouseholderQr().solve(Bz);

  vector<double> cx(6), cy(6), cz(6);
  for (int i = 0; i < 6; i++)
  {
    cx[i] = Cofx(i);
    cy[i] = Cofy(i);
    cz[i] = Cofz(i);
  }

  PolynomialTraj poly_traj;
  poly_traj.addSegment(cx, cy, cz, t);

  return poly_traj;
}
