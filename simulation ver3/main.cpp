#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <GL/glut.h>
#include <conio.h>
#include <algorithm>
#include <chrono>

#include "Difinition.h"
#include "Bitmap.h"
#include "Plain.h"
#include "Yuragi.h"
#include "Parameter.h"
#include "Pressure.h"

using namespace std;


//---------------------------------------------------------------------------------------------------
//　　clear
//　　Desc : 変数の初期化
//---------------------------------------------------------------------------------------------------
void clear(void) {

	for (int i = 0; i <= NN + 1; i++) {
		for (int j = 0; j <= NN + 1; j++) {
			for (int k = 0; k <= 1; k++) {

				c[i][j][k] = 0;
				water[i][j][k] = 0;
				dye[i][j][k] = 0;
				p[i][j][k] = P_AT;//大気圧1 atm
				weft[i][j][k] = 0; warp[i][j][k] = 0;
				Gap[i][j][k] = 0;
				for (int x = 0; x < 5; x++) AdjacentCellStatus[i][j][k][x] = 0;
				dCount[i][j][k] = 0;

			}
		}
	}

	cout << endl << "変数の初期化(clear)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　Cloth
//　　Desc : 布の初期値の設定
//---------------------------------------------------------------------------------------------------
void Cloth(void) {

	///////////////////////////////////////布の初期値の設定///////////////////////////////////////
	gapx = (25.4 - densityx * yarnx) / (densityx - 1);//横糸隙間の基準
	gapy = (25.4 - densityy * yarny) / (densityy - 1);//縦糸隙間の基準
	if (gapx > yarnx || gapy > yarny) cout << "＜注意＞糸の太さより隙間の方が大きい" << endl;
	if (gapx < 0 || gapy < 0) {
		cout << "隙間の大きさgapxまたはgapyが負の値のため...強制終了" << endl;
		exit(1);
	}
	//布内の各糸の本数計算，どっちもほぼかわらん
	//Nx = height_real / (yarnx + gapx); Ny = width_real / (yarny + gapy);//独自のちょいなぞな計算
	Nx = height_real * (densityx / 25.4);	Ny = width_real * (densityy / 25.4);//密度を用いて計算
	Nx = floor(Nx); Ny = floor(Ny);
	if ((2 * Nx + 1) > NN || (2 * Ny + 1) > NN) {
		cout << "布内の糸の本数が設定してある配列サイズをこえるため...強制終了" << endl;
		exit(1);
	}
	//描画用の変換（拡大）
	YARNx = converting_rate * yarnx; YARNy = converting_rate * yarny;
	GAPx = converting_rate * gapx; GAPy = converting_rate * gapy;
	width = converting_rate * width_real; height = converting_rate * height_real;
	//小数点以下切捨て
	width = floor(width); height = floor(height);
	//四捨五入の方法
	//Rx = floor(Rx + 0.5); Ry = floor(Ry + 0.5);

	//糸セルと隙間セルの一辺の長さを定義
	for (int i = 1; i <= (2 * Ny + 1); i++) {
		for (int j = 1; j <= (2 * Nx + 1); j++) {
			if (i % 2 == 0 && j % 2 == 1) {//横糸間の隙間
				gap[i][j][X] = yarny;//セル幅
				gap[i][j][Y] = gapx;//セル奥行，隙間の直径
				gap[i][j][Z] = (yarnx + yarny) / 2;//セル高さ
				yarn[i][j][X] = yarny;
				yarn[i][j][Y] = 0;
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 0) {//縦糸間の隙間
				gap[i][j][X] = gapy;//セル幅，隙間の半径
				gap[i][j][Y] = yarnx;//セル奥行
				gap[i][j][Z] = (yarnx + yarny) / 2;//セル高さ
				yarn[i][j][X] = 0;
				yarn[i][j][Y] = yarnx;
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 1) {//隙間
				gap[i][j][X] = gapy;//セル幅
				gap[i][j][Y] = gapx;//セル奥行
				gap[i][j][Z] = (yarnx + yarny) / 2;//セル高さ
				yarn[i][j][X] = 0;
				yarn[i][j][Y] = 0;
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 0 && j % 2 == 0) {//糸，布目
				gap[i][j][X] = 0;
				gap[i][j][Y] = 0;
				yarn[i][j][X] = yarny;//セル幅
				yarn[i][j][Y] = yarnx;//セル奥行
				yarn[i][j][Z] = (yarnx + yarny) / 2;//セル高さ
				capacity[i][j] = yarn[i][j][X] * yarn[i][j][Y] * yarn[i][j][Z];
			}
		}
	}

	cout << endl << "布の初期値設定(Cloth)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　Shibori
//　　Desc : 布の圧迫の設定，ゆくゆくは絞りの設定
//---------------------------------------------------------------------------------------------------
void Shibori(void) {

	//初期圧迫位置の設定
	double A1 = 0, A2 = 0, A3 = 0, A4 = 0;
	int n = 0;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			A1 = 0.9 * 2 * Ny - i; A2 = 1.1 * 2 * Ny - i;
			A3 = 0.9 * 2 * Ny - i; A4 = 1.1 * 2 * Ny - i;
			if (j > A1 && j < A2) {
				p[i][j][X] = p[i][j][Y] = P_MAX;
				p0_i[n] = i; p0_j[n] = j;//初期圧迫位置を記録
				n++;
			}
		}
	}

	//ガウシアン分布
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			//上で記録した圧迫位置のところのみ，圧力分散をガウス関数で計算
			for (int l = 0; l <= n; l++) {
				if (p0_i[l] == i && p0_j[l] == j) {
					Pressure(i, j);
				}
				p[i][j][Y] = p[i][j][X];//縦糸層と横糸層で同じ圧力分布
			}
		}
	}

	//圧力をもとに糸・隙間の高さを変化させる
	double ap;//圧力によるセルの膨張率，セルの幅・奥行きがどれくらい大きくなるか
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				if (p[i][j][k] > P_AT && p[i][j][k] < P_MAX) {
					gap[i][j][Z] = gap[i][j][Z] / p[i][j][k];
					yarn[i][j][Z] = yarn[i][j][Z] / p[i][j][k];
					ap = sqrt(p[i][j][k]);
					gap[i][j][X] = ap * gap[i][j][X]; gap[i][j][Y] = ap * gap[i][j][Y];
					yarn[i][j][X] = ap * yarn[i][j][X]; yarn[i][j][Y] = ap * yarn[i][j][Y];
				}
				else if (p[i][j][k] >= P_MAX) {
					p[i][j][k] = P_MAX;
					gap[i][j][Z] = 0;
					yarn[i][j][Z] = 0;
					ap = sqrt(p[i][j][k]);
					gap[i][j][X] = ap * gap[i][j][X]; gap[i][j][Y] = ap * gap[i][j][Y];
					yarn[i][j][X] = ap * yarn[i][j][X]; yarn[i][j][Y] = ap * yarn[i][j][Y];
				}
			}
		}
	}

	//ゆらぎと圧力に基づいてセルの許容量を再計算
	for (int i = 0; i <= (2 * Ny + 1); i++) {
		for (int j = 0; j <= (2 * Nx + 1); j++) {
			if (i % 2 == 0 && j % 2 == 1) {//横糸間の隙間
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 0) {//縦糸間の隙間
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 1) {//隙間
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 0 && j % 2 == 0) {//糸，布目
				capacity[i][j] = yarn[i][j][X] * yarn[i][j][Y] * yarn[i][j][Z];
			}
		}
	}

	cout << endl << "絞り設定(Shibori)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　Dye
//　　Desc : 染料の初期値の設定，染料の与える位置など
//---------------------------------------------------------------------------------------------------
void Dye(void) {
	//////////////////////////////////////染料の初期値の設定///////////////////////////////////////
	double rA;
	//染色前の布の状態を変更
	if (wet_or_dry == 1) {
		for (int i = 1; i <= 2 * Ny + 1; i++) {
			for (int j = 1; j <= 2 * Nx + 1; j++) {
				for (int k = 0; k <= 1; k++) {
					water[i][j][k] = capacity[i][j] * 1.0;
				}
			}
		}
	}
	//板で斜めに防染した位置よりも右上に円状に滴下
	double A;
	rA = 40;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				A = (i - 0.35 * 2 * Ny) * (i - 0.35 * 2 * Ny) + (j - 0.35 * 2 * Nx) * (j - 0.35 * 2 * Nx);
				if (p[i][j][k] >= P_AT && p[i][j][k] < P_MAX) {
					if (A <= rA * rA) {
						dye[i][j][k] += 0.1;
						water[i][j][k] += 0.1;
						c[i][j][k] = dye[i][j][k] / (water[i][j][k] + dye[i][j][k]);
						//cout << "(" << i << ", " << j << ")\t";
					}
				}
			}
		}
	}
	//中心に滴下
	/*rA = 60;
	cout << "滴下範囲の半径 rA = " << rA << endl;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				if (i % 2 == 1 || j % 2 == 1) {
					A = (i - Ny) * (i - Ny) + (j - Nx) * (j - Nx);
					if (A <= rA * rA) {
						dye[i][j][k] += 0.1;
						water[i][j][k] += 0.1;
						c[i][j][k] = dye[i][j][k] / (water[i][j][k] + dye[i][j][k]);
						dyewater += (dye[i][j][k] + water[i][j][k]);
					}
				}
			}
		}
	}*/
	cout << "染料溶液量：" << dyewater << endl;
	//中心1点のみに滴下
	/*for (int k = 0; k <= 1; k++) {
		dye[Ny + 1][Nx + 1][k] += 0.1;
		water[Ny + 1][Nx + 1][k] += 0.1;
		c[Ny + 1][Nx + 1][k] = dye[Ny + 1][Nx + 1][k] / (water[Ny + 1][Nx + 1][k] + dye[Ny + 1][Nx + 1][k]);
	}*/

	cout << "染料の初期値設定(Dye)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　ClothDraw
//　　Desc : 描画のときのセルの大きさとか布の大きさの設定
//---------------------------------------------------------------------------------------------------
void ClothDraw(void) {

	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			if (i % 2 == 0 && j % 2 == 0) {
				weft[i][j][X] = weft[i - 2][j][X] + 0.5 * gap[i - 3][j][X] + yarn[i - 2][j][X] + 0.5 * gap[i - 1][j][X];
				weft[i][j][Y] = weft[i][j - 2][Y] + yarn[i][j - 2][Y] + gap[i][j - 1][Y];
				warp[i][j][X] = warp[i - 2][j][X] + yarn[i - 2][j][X] + gap[i - 1][j][X];
				warp[i][j][Y] = warp[i][j - 2][Y] + 0.5 * gap[i][j - 3][Y] + yarn[i][j - 2][Y] + 0.5 * gap[i][j - 1][Y];
				//最左端と最下端の設定
				weft[2][j][X] = 0;
				weft[i][2][Y] = 0.5 * gap[1][1][Y];
				warp[2][j][X] = 0.5 * gap[1][1][X];
				warp[i][2][Y] = 0;
			}
			else if (i % 2 == 1 && j % 2 == 1) {
				Gap[i][j][X] = Gap[i - 2][j][X] + gap[i - 2][j][X] + gap[i - 1][j][X];
				Gap[i][j][Y] = Gap[i][j - 2][Y] + gap[i][j - 2][Y] + gap[i][j - 1][Y];
				//最左端と最下端の設定
				Gap[1][j][X] = -0.5 * gap[1][1][X];
				Gap[i][1][Y] = -0.5 * gap[1][1][Y];
			}
			/*weft[i][j][X] -= 0.1;
			weft[i][j][Y] -= 0.1;
			warp[i][j][X] -= 0.1;
			warp[i][j][Y] -= 0.1;
			Gap[i][j][X] -= 0.1;
			Gap[i][j][Y] -= 0.1;*/
		}
	}
	

	cout << "weft[2 * Ny][2 * Nx][X], [Y] = (" << weft[2 * Ny][2 * Nx][X] << ", " << weft[2 * Ny][2 * Nx][Y] << ")" << endl;
	cout << "warp[2 * Ny][2 * Nx][X], [Y] = (" << warp[2 * Ny][2 * Nx][X] << ", " << warp[2 * Ny][2 * Nx][Y] << ")" << endl;
	cout << "Gap[2 * Ny + 1][2 * Nx + 1][X], [Y] = (" << Gap[2 * Ny + 1][2 * Nx + 1][X] << ", " << Gap[2 * Ny + 1][2 * Nx + 1][Y] << ")" << endl;

	width = converting_rate * (Gap[2 * Ny + 1][2 * Nx + 1][X] + gap[2 * Ny + 1][2 * Nx + 1][X]);
	height = converting_rate * (Gap[2 * Ny + 1][2 * Nx + 1][Y] + gap[2 * Ny + 1][2 * Nx + 1][Y]);
	//width = floor(width); height = floor(height);

	cout << "描画の位置設定(ClothDraw)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　Yuragi
//　　Desc : ゆらぎの付与
//---------------------------------------------------------------------------------------------------
void Yuragi(void) {
	//ゆらぎの導入
	double rmin = NN;
	double mmax = 0;
	//ゆらぎの幅（yuragiwの最低値～最大値）を設定
	yuragixmax = 0;
	yuragiwmax = 0;
	for (int i = 1; i <= 2 * Ny + 1; i++) {//すべてのセルでゆらぎ計算，i = 0においてゆらぎは0
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			YuragiCal(i, j);
		}
	}
	yuragix_range_setting = 1.0;
	while (yuragixmax / yuragix_range_setting > yuragix_range) {
		yuragix_range_setting += 0.00001;
	}
	yuragiw_range_setting = 1.0;
	while (yuragiwmax / yuragiw_range_setting > yuragiw_range) {
		yuragiw_range_setting += 0.00001;
	}

	yuragixmax = 0, yuragixmin = NN;
	yuragiwmax = 0, yuragiwmin = NN;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			yuragix[i][j] = yuragix[i][j] / yuragix_range_setting;
			yuragixmax = (yuragix[i][j] >= yuragixmax) ? yuragix[i][j] : yuragixmax;
			yuragixmin = (yuragix[i][j] <= yuragixmin) ? yuragix[i][j] : yuragixmin;
			yuragiw[i][j] = yuragiw[i][j] / yuragiw_range_setting;
			yuragiwmax = (yuragiw[i][j] >= yuragiwmax) ? yuragiw[i][j] : yuragiwmax;
			yuragiwmin = (yuragiw[i][j] <= yuragiwmin) ? yuragiw[i][j] : yuragiwmin;
		}
	}

	for (int i = 1; i <= 2 * Ny + 1; i++) {//i = 0において計算するとゆらぎが0になる
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			yuragix[i][j] = yuragix[i][j] - (yuragixmax + yuragixmin) / 2;
			yuragiw[i][j] = yuragiw[i][j] - (yuragiwmax + yuragiwmin) / 2 + 1.0;
			weft[i][j][X] += yuragix[i][j];
			weft[i][j][Y] += yuragix[i][j];
			warp[i][j][X] += yuragix[i][j];
			warp[i][j][Y] += yuragix[i][j];
			Gap[i][j][X] += yuragix[i][j];
			Gap[i][j][Y] += yuragix[i][j];

			gap[i][j][X] *= yuragiw[i][j];
			gap[i][j][Y] *= yuragiw[i][j];
			yarn[i][j][X] *= yuragiw[i][j];
			yarn[i][j][Y] *= yuragiw[i][j];
		}
	}
	
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			weftsize[i][j][X] = 0.5 * gap[i - 1][j][X] + yarn[i][j][X] + 0.5 * gap[i + 1][j][X];
			weftsize[i][j][Y] = yarn[i][j][Y];
			warpsize[i][j][X] = yarn[i][j][X];
			warpsize[i][j][Y] = 0.5 * gap[i][j - 1][Y] + yarn[i][j][Y] + 0.5 * gap[i][j + 1][Y];
		}
	}
	double sukima;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			if (i % 2 == 0 && j % 2 == 0) {
				sukima = weft[i][j][X] - (weft[i - 2][j][X] + weftsize[i - 2][j][X]);
				if (sukima > 0) weftsize[i - 2][j][X] += sukima;
				sukima = warp[i][j][Y] - (warp[i][j - 2][Y] + warpsize[i][j - 2][Y]);
				if (sukima > 0) warpsize[i][j - 2][Y] += sukima;
			}
		}
	}
	
	cout << "ゆらぎ設定(Yuragi)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　Display1
//　　Desc : ウィンドウ1への描画
//---------------------------------------------------------------------------------------------------
void Display1(void) {
	auto startTime = std::chrono::system_clock::now();

	//画面クリア
	glClear(GL_COLOR_BUFFER_BIT);

	//染色計算
	double dyemax = 0;
	cout << endl << "染色計算中";
	cout << endl << "Press 'q' to get out of the loop" << endl << "Press 't' to get interim output" << endl;
	for (t = 1; t <= loop_times; t++) {
		dTerminator = 0;
		bcTerminator = 0;

		//着目セルの近傍のセルの状態(湿潤 or 乾燥)を確認
		for (int i = 1; i <= 2 * Ny + 1; i++) {
			for (int j = 1; j <= 2 * Nx + 1; j++) {
				for (int k = 0; k <= 1; k++) {
					//隣接セルの乾燥・湿潤を判別
					if (water[i + 1][j][k] >= capacity[i + 1][j] * 1.0) AdjacentCellStatus[i][j][k][0] = 1;
					if (water[i - 1][j][k] >= capacity[i - 1][j] * 1.0) AdjacentCellStatus[i][j][k][1] = 1;
					if (water[i][j + 1][k] >= capacity[i][j + 1] * 1.0) AdjacentCellStatus[i][j][k][2] = 1;
					if (water[i][j - 1][k] >= capacity[i][j - 1] * 1.0) AdjacentCellStatus[i][j][k][3] = 1;
					//隣接セルが圧迫されているか
					if (p[i + 1][j][k] >= P_MAX) AdjacentCellStatus[i][j][k][0] = 2;
					if (p[i - 1][j][k] >= P_MAX) AdjacentCellStatus[i][j][k][1] = 2;
					if (p[i][j + 1][k] >= P_MAX) AdjacentCellStatus[i][j][k][2] = 2;
					if (p[i][j - 1][k] >= P_MAX) AdjacentCellStatus[i][j][k][3] = 2;
					//他の層の状態は現在未考慮(2017/06/15)
					/*if (k == 0 && water[i][j][k + 1] >= capacity[i][j] * 1.0) AdjacentCellStatus[i][j][k][4] = 1;
					else if (k == 1 && water[i][j][k - 1] >= capacity[i][j] * 1.0) AdjacentCellStatus[i][j][k][4] = 1;*/
					dCount[i][j][k] = 0;
					bcCount[i][j][k] = 0;
					d_pequalCount[i][j][k] = 0;
				}
			}
		}

		if (d_pequalTerminator > 0) {//t = t-1のときの結果から
			for (int i = 1; i <= 2 * Ny + 1; i++) {
				for (int j = 1; j <= 2 * Nx + 1; j++) {
					for (int k = 0; k <= 1; k++) {
						if (p[i][j][k] > P_AT && p[i][j][k] < P_MAX) {
							PlateBoundaryCal(i, j, k);
						}
					}
				}
			}
			for (int i = 2 * Ny + 1; i >= 1; i--) {
				for (int j = 2 * Nx + 1; j >= 1; j--) {
					for (int k = 1; k >= 0; k--) {
						if (p[i][j][k] > P_AT && p[i][j][k] < P_MAX) {
							PlateBoundaryCal(i, j, k);
						}
					}
				}
			}
		}
		d_pequalTerminator = 0;

		for (int i = 1; i <= 2 * Ny + 1; i++) {
			for (int j = 1; j <= 2 * Nx + 1; j++) {
				for (int k = 0; k <= 1; k++) {
					PlainCal(i, j, k);
				}
			}
		}
		for (int i = 2 * Ny + 1; i >= 1; i--) {
			for (int j = 2 * Nx + 1; j >= 1; j--) {
				for (int k = 1; k >= 0; k--) {
					PlainCal(i, j, k);
				}
			}
		}

		for (int i = 1; i <= 2 * Ny + 1; i++) {
			for (int j = 1; j <= 2 * Nx + 1; j++) {
				for (int k = 0; k <= 1; k++) {
					if (dCount[i][j][k] != 0) dTerminator++;//すべてのセルで拡散計算されなければdTerminator = 0
					if (bcCount[i][j][k] != 0) bcTerminator++;//すべてのセルでBurasの式・毛細管作用の計算がされなければbcTerminator = 0
					if (d_pequalCount[i][j][k] != 0) d_pequalTerminator++;
				}
			}
		}
		if (t % 100 == 0) {
			auto endTime = std::chrono::system_clock::now();
			auto timeSpan = endTime - startTime;
			cout << "t = " << t << " (" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms])" << 
				endl << "　dTerminator = " << dTerminator << ", bcTerminator = " << bcTerminator << ", d_pequalTerminator = " << d_pequalTerminator << endl;
		}
		if (dTerminator == 0 && bcTerminator == 0 && d_pequalTerminator == 0) {
			cout << "t = " << t << endl;
			break;
		}

		//特定のキーが押下されたときの処理
		char buf;
		if (_kbhit()) {
			buf = _getch();
			//qキーが押下されたときにループ抜け出し
			if (buf == 'q') {
				cout << "t = " << t << "で強制終了" << endl;
				break;
			}
			//tキーが押下されたときに途中結果画像出力
			if (buf == 't') {
				for (int i = 1; i <= 2 * Ny + 1; i++) {
					for (int j = 1; j <= 2 * Nx + 1; j++) {
						for (int k = 0; k <= 1; k++) {
							if (capacity[i][j] != 0) {
								dyeDraw[i][j][k] = dye[i][j][k] / capacity[i][j] * 2;
							}
						}
					}
				}
				for (int i = 1; i <= 2 * Ny + 1; i++) {
					for (int j = 1; j <= 2 * Nx + 1; j++) {
						for (int k = 0; k <= 1; k++) {
							DrawGapPlain(i, j, k);
						}
					}
				}
				for (int i = 1; i <= 2 * Ny + 1; i++) {
					for (int j = 1; j <= 2 * Nx + 1; j++) {
						for (int k = 0; k <= 1; k++) {
							DrawYarnPlain(i, j, k);
						}
					}
				}
				WriteBitmap("Simulation output(途中経過).bmp", width, height);
				auto endTime = std::chrono::system_clock::now();
				auto timeSpan = endTime - startTime;
				cout << "途中経過画像出力: t = " << t << " (" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms])" << endl;
			}
		}

	}

	cout << "染色...COMPLETE" << endl;

	//染料量計算
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				if (capacity[i][j] != 0) {
					dyeDraw[i][j][k] = dye[i][j][k] / capacity[i][j] * 2;
				}
			}
		}
	}

	//描画
	cout << endl << "描画" << endl;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				DrawGapPlain(i, j, k);
			}
		}
	}
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				DrawYarnPlain(i, j, k);
			}
		}
	}

	cout << "描画...COMPLETE" << endl;
	WriteBitmap("Simulation output.bmp", width, height);
	cout << "保存...COMPLETE" << endl;

	glFlush();

	auto endTime = std::chrono::system_clock::now();
	auto timeSpan = endTime - startTime;
	cout << "処理時間:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　Display2
//　　Desc : ウィンドウ2への描画
//---------------------------------------------------------------------------------------------------
void Display2(void) {

	//画面クリア
	glClear(GL_COLOR_BUFFER_BIT);

	//圧力分布
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
			}
		}
	}
	glFlush();
}

//---------------------------------------------------------------------------------------------------
//　　resize
//　　Desc : 
//---------------------------------------------------------------------------------------------------
void resize(int w, int h) {
	/* ウィンドウ全体をビューポートにする */
	glViewport(0, 0, w, h);
	/* 変換行列の初期化 */
	glLoadIdentity();
	//表示する座標範囲を(0, 0, -1)~(400, 400, 1)とする
	glOrtho(0, w, 0, h, -1.0, 1.0);
}

//----------------------------------------------------------------------------------------------------
//　　Initialize
//　　Desc : 初期化処理，背景色
//----------------------------------------------------------------------------------------------------
void Initialize(void) {
	//glClearColor(R, G, B, α値：不透明度（0 で透明, 1 で不透明）);
	glClearColor(1.0, 1.0, 1.0, 1.0);
}

//----------------------------------------------------------------------------------------------------
//　　main
//　　Desc : メインエントリーポイント
//----------------------------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
	
	int WinID[2];

	Parameter();//入力パラメータ
	clear();
	Cloth();
	ClothDraw();//ゆらぎなしで描画位置を計算
	Shibori();
	Yuragi();	
	Dye();

	cout << "隙間：" << endl << "　縦糸間 gapy = " << gapy << endl;
	cout << "　横糸間 gapx = " << gapx << endl;
	cout << "布内の糸の本数：" << endl << "　縦糸 Ny = " << Ny << "　(セル数 2Ny+1 = " << 2 * Ny + 1 << ")"
		<< endl << "　横糸 Nx = " << Nx << "　(セル数 2Nx+1 = " << 2 * Nx + 1 << ")" << endl;
	cout << "セルの標準体積：" << endl << "　糸：" << yarnx * yarny * (yarnx + yarny) / 2 << endl;
	cout << "　隙間(縦横)：" << gapx * yarny * (yarnx + yarny) / 2 << endl;
	cout << "布の許容量(体積，ゆらぎなし)：" << ((yarny + gapy) * Nx) * ((yarnx + gapx) * Ny) * (yarnx + yarny) << endl;
	/*cout << endl << "YARNx = " << YARNx << ", YARNy = " << YARNy << ", GAPx = " << GAPx << ", GAPy = " << GAPy << endl;*/
	cout << "width = " << width << ", height = " << height << endl;
	cout << "loop times = ";	cout << loop_times << endl;
	//cin >> loop_times;

	glutInit(&argc, argv);//glutの初期化
	//ウィンドウ1
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);//値はピクセル
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	WinID[0] = glutCreateWindow("Tie Dye simulation");
	glutDisplayFunc(Display1);
	glutReshapeFunc(resize);
	glutMouseFunc(mouse);
	Initialize();

	//ウィンドウ2
	/*glutInitWindowPosition(1000, 0);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	WinID[1] = glutCreateWindow("Simulation Check");
	glutDisplayFunc(Display2);
	glutReshapeFunc(resize);
	glutMouseFunc(mouse);
	Initialize();*/

	glutMainLoop();
	return 0;

}