#pragma once
#define NN 1000//布の一辺の長さ（糸+隙間）の最大値，2*Nxや2*Nyより大きくなくてはならない
int Nx, Ny;//布内の横糸の本数，縦糸の本数
int width, height;//ウィンドウ大きさ・座標
int converting_rate;//シミュレーション（計算，実際の布サイズ）と画面表示の変換比率
int width_real, height_real;//実際の布のサイズ，10 cmとか

double c[NN + 2][NN + 2][2];//セルの染料の濃度（0：緯糸, 1：経糸）（t)
//濃度の拡散計算のとき使用
double dcN;//糸 → 隙間
double dcS;
double dcW;
double dcE;
double dcN_gap;//隙間 → 隙間
double dcS_gap;
double dcW_gap;
double dcE_gap;
double dcN_yarn;//隙間 → 糸
double dcS_yarn;
double dcW_yarn;
double dcE_yarn;

double dye[NN + 2][NN + 2][2];//染料量
double ddye;//染料の移動量
double ddye1;//染料の移動量

double water[NN + 2][NN + 2][2];//水分量
double dwater;//移動する水分，セルの飽和量を超えた分
double dwater1;//dwaterを1セル分に調整，Plain.h参照

double dyewater;//染料溶液量
double ddyewater;//浸透計算で移動する染料溶液量
double ddyewater1;//浸透計算で移動する染料溶液量

//Burasの式
double dq;//時間当たりの吸水量
double dq1;//時間当たりの吸水量
double initialVelocity = 0.1;//初期吸水速度
double elapsedTime;//経過時間

//拡散係数
double D;
//隙間 → 糸への拡散
double DyN;//着目セルの↑
double DyS;//↓
double DyW;//←
double DyE;//→
double Dy0;//乾燥糸へ
double Dy1;//湿潤糸へ
//隙間 → 隙間 or 糸への拡散
double DgN;//↑
double DgS;//↓
double DgW;//←
double DgE;//→
double Dg0;//隙間 → 隙間，乾燥隙間へ
double Dg1;//隙間 → 隙間，湿潤隙間へ
double Dg2;//糸 → 隙間，湿潤隙間へ

//隙間 → 糸における水分浸透の係数
double Dw;

//ゆらぎ
double yuragix[NN][NN];
double yuragiy[NN][NN];
double yuragiw[NN][NN];
double ax, bx, cx;
double ay, by, cy;
double aw, bw, cw;
int tajudox, tajudow;
double nx, ex, nw, ew;
double fx1, fx2, fy1, fy2, fw1, fw2;
double yuragixmax, yuragixmin;
double yuragiymax, yuragiymin;
double yuragiwmax, yuragiwmin;
double yuragix_range;//ゆらぎの幅，1-yuragi_range ~ 1+yuragi_range，0.8~1.2みたいな
double yuragiw_range;//ゆらぎの幅，1-yuragi_range ~ 1+yuragi_range，0.8~1.2みたいな
double yuragix_range_setting;//ゆらぎ量をyuragi_rangeにするため
double yuragiw_range_setting;//ゆらぎ量をyuragi_rangeにするため

//圧迫
double p[NN + 2][NN + 2][2];//圧力
int p0_i[NN * 100], p0_j[NN * 100];//圧迫位置，最初に布を圧迫しているところ

//毛細管作用
int Temperature;//温度
double SurfaceTension;//表面張力[N/m]
const double ContactAngle = 0;//接触角θ
const double LiquidDensity = 1000;//液体の密度[kg/m^3]
double Viscosity;//液体の粘度（粘性係数）[Pa*s]
double r;//毛細管作用の浸透距離を計算する際の隙間の半径
double rx;//x(i)方向の浸透距離計算のときの隙間半径
double ry;//y(j)方向の浸透距離計算のときの隙間半径
double h;//毛細管作用による浸透距離
double h1;//計算用
double hx;//毛細管作用によるx(i)方向の浸透距離
double hx1;//計算用
double hy;//毛細管作用によるy(j)方向の浸透距離
double hy1;//計算用
int hpositivex, hnegativex, hpositivey, hnegativey;//浸透セル数をカウント


//布の設定
double densityx;//横糸の密度（130~253，1インチ=25.4cmより）
double densityy;//縦糸の密度
double gap[NN + 2][NN + 2][3];//糸間の隙間の大きさ[mm]，0：セル幅x，1：セル奥行y，2：セル高さz
double gapx;//横糸間の隙間の直径の基準値[mm]
double gapy;//縦糸間の隙間の直径の基準値[mm]
double yarn[NN + 2][NN + 2][3];//糸の太さ = セルの一辺の長さ[mm]，0：セル幅x，1：セル奥行y，2：セル高さz
double yarnx;//横糸の太さの基準値（いわゆる番手より得られる太さ）[mm]
double yarny;//縦糸の太さの基準値（いわゆる番手より得られる太さ）[mm]
double capacity[NN + 2][NN + 2];//セルの飽和量（許容量，体積）[mm^3]


//計算用の変数
double dt = 0.005;//1タイムステップに経過する時間0.005
double c0 = 0;//初期染料濃度,0.4
double dye0 = 0.2;//初期染料量,0.2
double water0 = 0.3;//初期水分量,0.3
int t = 0;//計算を回した回数，時間
int loop_times = 0;
const int X = 0;
const int Y = 1;
const int Z = 2;
int WetDryFlag[NN + 2][NN + 2][2][5];//値が0のとき乾燥・1のとき湿潤，0:着目セルの右のセル，1:左，2:上，3:下，4:他の層(隙間において)
int dTerminator = 0;//染色終了条件(拡散計算)
int bcTerminator = 0;//染色終了条件(Burasの式と毛細管作用)
int dCount[NN + 2][NN + 2][2];//拡散計算が起きているか確認用
int bcCount[NN + 2][NN + 2][2];//Burasの式と毛細管作用の計算が起きているか確認用
int wet_or_dry;//初期状態において布は湿潤か乾燥か(0：乾燥，1：湿潤)


//描画用
int YARNx, YARNy;//【描画用】糸セルの一辺（直径）
double GAPx, GAPy;//【描画用】隙間セルの一辺（直径）
int I, J;
double weft[NN + 2][NN + 2][2];//描画開始座標，横糸セルの左下角の座標値，0：X座標，1：Y座標
double warp[NN + 2][NN + 2][2];//描画開始座標，縦糸セルの左下角の座標値，0：X座標，1：Y座標
double Gap[NN + 2][NN + 2][2];//描画開始座標，布の隙間部分の四角形の左下角の座標値，0：X座標，1：Y座標
double weft1[NN + 2][NN + 2][2];//描画開始座標，横糸セルの左下角の座標値，0：X座標，1：Y座標
double warp1[NN + 2][NN + 2][2];//描画開始座標，縦糸セルの左下角の座標値，0：X座標，1：Y座標
double Gap1[NN + 2][NN + 2][2];//描画開始座標，布の隙間部分の四角形の左下角の座標値，0：X座標，1：Y座標
double dyeDraw[NN + 2][NN + 2][2];//描画用染料値

//織り方
int weavepattern = 0;//Plain = 1, Twill = 2，Sateen = 3

//---------------------------------------------------------------------------------------------------
//　　InitRand
//　　Desc : 乱数の初期化
//---------------------------------------------------------------------------------------------------
inline void InitRand(void)
{
	srand((unsigned int)time(NULL));
}

//---------------------------------------------------------------------------------------------------
//　　Rand_mod
//　　Desc : 乱数，使いやすいように調整，a:乱数範囲（割る数，あまり），b:乱数の桁調整
//---------------------------------------------------------------------------------------------------
inline int Rand_mod(int a, double b)
{
	return rand() % a * b;//rand()は0~RAND_MAXまで
}

//---------------------------------------------------------------------------------------------------
//　　Rand_range
//　　Desc : 乱数，乱数範囲：min~max
//---------------------------------------------------------------------------------------------------
inline int Rand_range(double min, double max)
{
	return min + (double)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

