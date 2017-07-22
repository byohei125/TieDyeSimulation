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
				p[i][j][k] = 1;//大気圧1 atm
				weft[i][j][k] = 0; warp[i][j][k] = 0;
				Gap[i][j][k] = 0;
				for (int x = 0; x < 5; x++) WetDryFlag[i][j][k][x] = 0;
				dCount[i][j][k] = 0;

			}
		}
	}

	cout << endl << "変数の初期化...COMPLETE" << endl;

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
	for (int i = 0; i <= (2 * Ny + 1); i++) {
		for (int j = 0; j <= (2 * Nx + 1); j++) {
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
	

	cout << endl << "布の初期値設定...COMPLETE" << endl;

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
				p[i][j][X] = 1.5; p[i][j][Y] = 1.5;
				p0_i[n] = i; p0_j[n] = j;//初期圧迫位置を記録
				n++;
			}
			/*else if (j < A3 && j > A4) {
				p[i][j][X] = 1.3; p[i][j][Y] = 1.3;
				p0_i[n] = i; p0_j[n] = j;
				n++;
			}*/
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
				if (p[i][j][k] > 1.0 && p[i][j][k] < 1.5) {
					gap[i][j][Z] = gap[i][j][Z] / p[i][j][k];
					yarn[i][j][Z] = yarn[i][j][Z] / p[i][j][k];
					ap = sqrt(p[i][j][k]);
					gap[i][j][X] = ap * gap[i][j][X]; gap[i][j][Y] = ap * gap[i][j][Y];
					yarn[i][j][X] = ap * yarn[i][j][X]; yarn[i][j][Y] = ap * yarn[i][j][Y];
				}
				else if (p[i][j][k] >= 1.5) {
					p[i][j][k] = 1.5;
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
	/*for (int i = 0; i <= (2 * Ny + 1); i++) {
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
	}*/

}

//---------------------------------------------------------------------------------------------------
//　　Dye
//　　Desc : 染料の初期値の設定，染料の与える位置など
//---------------------------------------------------------------------------------------------------
void Dye(void) {
	//////////////////////////////////////染料の初期値の設定///////////////////////////////////////
	cout << endl << "染料の初期値設定" << endl;

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
	/*for (int i = Ny; i <= 2 * Ny + 1; i++) {
		for (int j = Nx; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				A = (i - 0.65 * 2 * Ny) * (i - 0.65 * 2 * Ny) + (j - 0.65 * 2 * Nx) * (j - 0.65 * 2 * Nx);
				if (A <= rA * rA) {
					dye[i][j][k] += 0.1;
					water[i][j][k] += 0.1;
					c[i][j][k] = dye[i][j][k] / (water[i][j][k] + dye[i][j][k]);
				}
			}
		}
	}*/
	//中心に滴下
	//rA = floor((Nx + Ny)/2 /2);
	rA = 60;//80
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
						//cout << "(" << i << ", " << j << ")\t";
					}
				}

			}
		}
	}
	cout << "染料溶液量：" << dyewater << endl;
	//中心1点のみに滴下
	/*for (int k = 0; k <= 1; k++) {
		dye[Ny + 1][Nx + 1][k] += 0.1;
		water[Ny + 1][Nx + 1][k] += 0.1;
		c[Ny + 1][Nx + 1][k] = dye[Ny + 1][Nx + 1][k] / (water[Ny + 1][Nx + 1][k] + dye[Ny + 1][Nx + 1][k]);
	}*/
	
	cout << "染料の初期値設定...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//　　ClothDraw
//　　Desc : 描画のときのセルの大きさとか布の大きさの設定
//---------------------------------------------------------------------------------------------------
void ClothDraw(void) {

	//初期値
	weft[2][2][X] = 0; weft[2][2][Y] = 0.5 * gap[1][1][Y];
	warp[2][2][X] = 0.5 * gap[1][1][X]; warp[2][2][Y] = 0;
	Gap[1][1][X] = -0.5 * gap[1][1][X]; Gap[1][1][Y] = -0.5 * gap[1][1][Y];

	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			if (i % 2 == 0 && j % 2 == 0) {
				if (i > 3 || j > 3) {
					weft[i][j][X] = weft[i - 2][j][X] + 0.5 * gap[i - 3][j][X] + yarn[i - 2][j][X] + 0.5 * gap[i - 1][j][X];
					weft[i][j][Y] = weft[i][j - 2][Y] + yarn[i][j - 2][Y] + gap[i][j - 1][Y];
					warp[i][j][X] = warp[i - 2][j][X] + yarn[i - 2][j][X] + gap[i - 1][j][X];
					warp[i][j][Y] = warp[i][j - 2][Y] + 0.5 * gap[i][j - 3][Y] + yarn[i][j - 2][Y] + 0.5 * gap[i][j - 1][Y];
				}
			}
			else if (i % 2 == 1 && j % 2 == 1) {
				Gap[i][j][X] = Gap[i - 2][j][X] + gap[i - 2][j][X] + gap[i - 1][j][X];
				Gap[i][j][Y] = Gap[i][j - 2][Y] + gap[i][j - 2][Y] + gap[i][j - 1][Y];
			}
		}
	}

	cout << "weft[2 * Ny][2 * Nx][X], [Y] = (" << weft[2 * Ny][2 * Nx][X] << ", " << weft[2 * Ny][2 * Nx][Y] << ")" << endl;
	cout << "warp[2 * Ny][2 * Nx][X], [Y] = (" << warp[2 * Ny][2 * Nx][X] << ", " << warp[2 * Ny][2 * Nx][Y] << ")" << endl;
	cout << "Gap[2 * Ny + 1][2 * Nx + 1][X], [Y] = (" << Gap[2 * Ny + 1][2 * Nx + 1][X] << ", " << Gap[2 * Ny + 1][2 * Nx + 1][Y] << ")" << endl;

	width = converting_rate * (Gap[2 * Ny + 1][2 * Nx + 1][X] + gap[2 * Ny + 1][2 * Nx + 1][X]);
	height = converting_rate * (Gap[2 * Ny + 1][2 * Nx + 1][Y] + gap[2 * Ny + 1][2 * Nx + 1][Y]);
	//width = floor(width); height = floor(height);

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
	yuragiwmax = 0, yuragiwmin = NN;
	int len = 0;
	for (int i = 1; i <= 2 * Ny + 1; i++) {//すべてのセルでゆらぎ計算，i = 0においてゆらぎは0
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			YuragiCal(i, j);
			yuragi_sort[len] = yuragiw[i][j];
			len++;
		}
	}

	//ゆらぎ量の中央値を求める
	sort(yuragi_sort, yuragi_sort + len);
	if (len % 2 == 0) yuragi_med = (yuragi_sort[len / 2 - 1] + yuragi_sort[len / 2]) / 2;
	else yuragi_med = yuragi_sort[len / 2];
	cout << "yuragi_med = " << yuragi_med << endl;

	yuragi_range_setting = 1.0;
	while (yuragiwmax / yuragi_range_setting > yuragi_range) {
		yuragi_range_setting += 0.00001;
	}

	for (int i = 1; i <= 2 * Ny + 1; i++) {//i = 0において計算するとゆらぎが0になる
		for (int j = 1; j <= 2 * Nx + 1; j++) {

			//if (i % 2 == 0 && j % 2 == 0) {//糸，布目
			//	yuragiw[i][j] = yuragiw[i][j] / yuragi_range_setting + 1.0;
			//	gap[i - 1][j][X] = gap[i - 1][j][X] + yarn[i][j][X] - yarny * yuragiw[i][j];
			//	if (gap[i - 1][j][X] < 0) {
			//		gap[i - 1][j][X] = 0;
			//	}
			//	yarn[i][j][X] = yarny * yuragiw[i][j];
			//	gap[i][j - 1][Y] = gap[i][j - 1][Y] + yarn[i][j][Y] - yarnx * yuragiw[i][j];
			//	if (gap[i][j - 1][Y] < 0) {
			//		gap[i][j - 1][Y] = 0;
			//	}
			//	yarn[i][j][Y] = yarnx * yuragiw[i][j];
			//}
			yuragiw[i][j] = yuragiw[i][j] / yuragi_range_setting + 1.0 - yuragi_med;
			gap[i][j][X] *= yuragiw[i][j];
			gap[i][j][Y] *= yuragiw[i][j];
			yarn[i][j][X] *= yuragiw[i][j];
			yarn[i][j][Y] *= yuragiw[i][j];

			/*if (i % 2 == 0 && j % 2 == 0) {
				if (i > 3 || j > 3) {
					weft[i][j][X] *= yuragiw[i][j];
					weft[i][j][Y] *= yuragiw[i][j];
					warp[i][j][X] *= yuragiw[i][j];
					warp[i][j][Y] *= yuragiw[i][j];
				}
			}
			else if (i % 2 == 1 && j % 2 == 1) {
				Gap[i][j][X] *= yuragiw[i][j];
				Gap[i][j][Y] *= yuragiw[i][j];
			}*/

		}
	}

	//横方向の糸の太さのつながりを確認
	/*for (int j = 1; j <= 2 * Nx + 1; j++) {
		for (int i = 1; i <= 2 * Ny + 1; i++) {	
			if (i % 2 == 0 && j % 2 == 0) {
				cout << "(" << i - 1 << "," << j << "), (" << i << "," << j << ")：" << gap[i - 1][j][Y] << ", " << yarn[i][j][Y] << endl;
			}
		}
	}*/

}

//---------------------------------------------------------------------------------------------------
//　　Display1
//　　Desc : ウィンドウ1への描画
//---------------------------------------------------------------------------------------------------
void Display1(void) {
	
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
					if (water[i + 1][j][k] >= capacity[i + 1][j] * 1.0) WetDryFlag[i][j][k][0] = 1;
					if (water[i - 1][j][k] >= capacity[i - 1][j] * 1.0) WetDryFlag[i][j][k][1] = 1;
					if (water[i][j + 1][k] >= capacity[i][j + 1] * 1.0) WetDryFlag[i][j][k][2] = 1;
					if (water[i][j - 1][k] >= capacity[i][j - 1] * 1.0) WetDryFlag[i][j][k][3] = 1;
					//他の層の状態は現在未考慮(2017/06/15)
					/*if (k == 0 && water[i][j][k + 1] >= capacity[i][j] * 1.0) WetDryFlag[i][j][k][4] = 1;
					else if (k == 1 && water[i][j][k - 1] >= capacity[i][j] * 1.0) WetDryFlag[i][j][k][4] = 1;*/
					dCount[i][j][k] = 0;
					bcCount[i][j][k] = 0;
				}
			}
		}

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
		
		if (t == (loop_times * 0.3) || t == (loop_times * 0.6) || t == (loop_times * 0.9)) cout << ".";

		for (int i = 1; i <= 2 * Ny + 1; i++) {
			for (int j = 1; j <= 2 * Nx + 1; j++) {
				for (int k = 0; k <= 1; k++) {
					if (dCount[i][j][k] != 0) dTerminator++;//すべてのセルで拡散計算されなければdTerminator = 0
					if (bcCount[i][j][k] != 0) bcTerminator++;//すべてのセルでBurasの式・毛細管作用の計算がされなければbcTerminator = 0
					/*if (isnan(c[i][j][k]) == 1) cout << "t = " << t << ", c[" << i << "][" << j << "][" << k << "] = " << c[i][j][k] << endl;
					if (isnan(water[i][j][k]) == 1) cout << "t = " << t << ", water[" << i << "][" << j << "][" << k << "] = " << water[i][j][k] << endl;
					if (isnan(dye[i][j][k]) == 1) cout << "t = " << t << ", dye[" << i << "][" << j << "][" << k << "] = " << dye[i][j][k] << endl;*/
				}
			}
		}
		if (t % 100 == 0) cout << "dTerminator = " << dTerminator << ", bcTerminator = " << bcTerminator << " (t = " << t << ")" << endl;
		if (dTerminator == 0 && bcTerminator == 0) {
			cout << "t = " << t << endl;
			break;
		}

		//qキーが押下されたときにループ抜け出し
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
							dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
							if (i % 2 == 0 && j % 2 == 0) {
								dyemax = (dye[i][j][k] >= dyemax) ? dye[i][j][k] : dyemax;
							}
						}
					}
				}
				for (int i = 1; i <= 2 * Ny + 1; i++) {
					for (int j = 1; j <= 2 * Nx + 1; j++) {
						for (int k = 0; k <= 1; k++) {
							dyeDraw[i][j][k] = dye[i][j][k] / dyemax;
						}
					}
				}
				for (int i = 1; i <= 2 * Ny + 1; i++) {
					for (int j = 1; j <= 2 * Nx + 1; j++) {
						for (int k = 0; k <= 1; k++) {
							DrawPlain(i, j, k);
						}
					}
				}
				WriteBitmap("Simulation output(途中).bmp", width, height);
			}
		}

	}

	cout << "染色...COMPLETE" << endl;
	
	//染料量計算
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				if (i % 2 == 0 && j % 2 == 0) {
					dyemax = (dye[i][j][k] >= dyemax) ? dye[i][j][k] : dyemax;
				}
			}
		}
	}
	cout << "dyemax = " << dyemax;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				//dye[i][j][k] /= dyemax;
				dyeDraw[i][j][k] = dye[i][j][k] / dyemax;
			}
		}
	}

	//描画
	cout << endl << "描画" << endl;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				DrawPlain(i, j, k);//セル大きさに沿って、セルを配置
			}
		}
	}
	
	cout << "描画...COMPLETE" << endl;
	WriteBitmap("Simulation output.bmp", width, height);
	cout << "保存...COMPLETE" << endl;

	glFlush();

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
				PressureCheck(i, j, k);
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
	cout << "アンケート用？（0：Yes，1：No）"; cin >> parameter_setting;
	Parameter(parameter_setting);//入力パラメータ
	
	clear();
	Cloth();
	Yuragi();
	ClothDraw();
	//Shibori();
	Dye();
	//Yuragi();

	cout << "隙間：" << endl << "　縦糸間 gapy = " << gapy << endl << "　　　　 gapy = " << gapy << endl;
	cout << "　縦糸間 gapx = " << gapx << endl << "　　　　 gapx = " << gapx << endl;
	cout << "布内の糸の本数：" << endl << "　縦糸 Ny = " << Ny  << "　(セル数 2Ny+1 = " << 2 * Ny + 1 << ")"
		<< endl << "　横糸 Nx = " << Nx << "　(セル数 2Nx+1 = " << 2 * Nx + 1 << ")" << endl;
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