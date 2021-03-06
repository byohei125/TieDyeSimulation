#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <Windows.h>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <GL/glut.h>
#include "Difinition.h"
using namespace std;

//---------------------------------------------------------------------------------------------------
//　　PlateBoundaryCal
//　　Desc :	板で防染されたときの，板の境界部における染料の流れ・浸透
//---------------------------------------------------------------------------------------------------
void PlateBoundaryCal(int i, int j, int k) {

	//着目セルに対して，等圧力のセルを抽出
	const int nn = 9;
	int pequal_i[nn] = {}, pequal_j[nn] = {};
	int n = 0;
	for (int l = -1; l <= 1; l++) {
		for (int m = -1; m <= 1; m++) {
			if (p[i][j][k] == p[i + l][j + m][k]) {
				pequal_i[n] = i + l;
				pequal_j[n] = j + m;
				n++;
			}
		}
	}
	if (n > 0 && c[i][j][k] > 0) {
		double dc_pequal[nn] = {};
		for (int l = 0; l < n; l++) {
			//着目セルと等圧力のセルとの距離の算出
			double d_pequal = 0;
			if ((i % 2 == 0 && j % 2 == 1) || (i % 2 == 1 && j % 2 == 0)) {//縦横隙間
				if (i != pequal_i[n] && j != pequal_j[n]) {//斜め方向に等圧力セルがあるとき
					d_pequal = sqrt(pow((gap[i][j][X] / 2), 2)  + pow((gap[i][j][Y] / 2), 2)) + sqrt(pow((gap[pequal_i[l]][pequal_j[l]][X] / 2), 2) * pow((gap[pequal_i[l]][pequal_j[l]][Y] / 2), 2));
				}
				else if (i != pequal_i[n]) {//横方向
					d_pequal = (gap[i][j][X] + gap[pequal_i[l]][pequal_j[l]][X]) / 2;
				}
				else if (j != pequal_j[n]) {//縦方向
					d_pequal = (gap[i][j][Y] + yarn[pequal_i[l]][pequal_j[l]][Y]) / 2;
				}
			}
			else if (i % 2 == 1 && j % 2 == 1) {//隙間
				if (i != pequal_i[n] && j != pequal_j[n]) {//斜め方向に等圧力セルがあるとき
					d_pequal = sqrt((gap[i][j][X] / 2) * (gap[i][j][Y] / 2)) + sqrt((yarn[pequal_i[l]][pequal_j[l]][X] / 2) * (yarn[pequal_i[l]][pequal_j[l]][Y] / 2));
				}
				else if (i != pequal_i[n]) {//横方向
					d_pequal = (gap[i][j][X] + gap[pequal_i[l]][pequal_j[l]][X]) / 2;
				}
				else if (j != pequal_j[n]) {//縦方向
					d_pequal = (gap[i][j][Y] + gap[pequal_i[l]][pequal_j[l]][Y]) / 2;
				}
			}
			else {//糸，糸から染み出すか？
				break;
			}

			//double dX1, dY1, dX2, dY2;
			//if (i % 2 == 0 && j % 2 == 0) {//糸
			//	dX1 = (weft[i][j][X] + warp[i][j][X]) / 2; dY1 = (weft[i][j][Y] + warp[i][j][Y]) / 2;				
			//}
			//else {//隙間
			//	dX1 = Gap[i][j][X]; dY1 = Gap[i][j][Y];
			//}
			//if (pequal_i[l] % 2 == 0 && pequal_j[l] % 2 == 0) {//糸
			//	dX2 = (weft[pequal_i[l]][pequal_j[l]][X] + warp[pequal_i[l]][pequal_j[l]][X]) / 2; dY2 = (weft[pequal_i[l]][pequal_j[l]][Y] + warp[pequal_i[l]][pequal_j[l]][Y]) / 2;
			//}
			//else {//隙間
			//	dX2 = Gap[pequal_i[l]][pequal_j[l]][X]; dY2 = Gap[pequal_i[l]][pequal_j[l]][Y];
			//}
			//d_pequal = sqrt(pow((dX2 - dX1), 2) + pow((dY2 - dY1), 2));

			dc_pequal[l] = dt *  D_pequal * ((c[i][j][k] - c[pequal_i[l]][pequal_j[l]][k]) / (pow(d_pequal, 2)));
			//cout << "dc_pequal[l] = " << dc_pequal[l] << endl;
			if (dc_pequal[l] > 0) {
				c[pequal_i[l]][pequal_j[l]][k] += dc_pequal[l];
				dye[pequal_i[l]][pequal_j[l]][k] = (c[pequal_i[l]][pequal_j[l]][k] / (1 - c[pequal_i[l]][pequal_j[l]][k])) * water[pequal_i[l]][pequal_j[l]][k];
				c[i][j][k] -= dc_pequal[l];
				dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
			}

			//染色終了条件
			double dThreshold = pow(10, -3) * 3.5;//pow(10, -4) * 3.5
			if (dc_pequal[l] > dThreshold) d_pequalCount[i][j][k]++;
		}
	}

}


//---------------------------------------------------------------------------------------------------
//　　PlainCal
//　　Desc : 平織，染色の計算の判断，毛細管なのか拡散なのか
//---------------------------------------------------------------------------------------------------
void PlainCal(int i, int j, int k) {

	hx = hx1 = hy = hy1 = 0;//値の初期化
	//セルの水分量が飽和量を超えた時点で（どんなに少量でも）隣のセルに水分を受け渡す
	hpositivex = hnegativex = hpositivey = hnegativey = 0;
	r = 0;
	dq = dq1 = 0;

	//隙間 → 隙間，隙間 → 糸における拡散移動量を計算
	dcN = dcS = dcW = dcE = 0;
	dcN_yarn = dcS_yarn = dcW_yarn = dcE_yarn = 0;
	dcN_gap = dcS_gap = dcW_gap = dcE_gap = 0;
	DyN = DyS = DyW = DyE = Dy1;
	DgN = DgS = DgW = DgE = Dg1;
	//x(i)方向の毛細管作用による浸透距離
	rx = (gap[i][j][Y] <= gap[i][j][Z]) ? (gap[i][j][Y] / 2) : (gap[i][j][Z] / 2);
	hx = sqrt(((rx * 0.001) * SurfaceTension * cos(ContactAngle) * dt) / (2 * Viscosity)) * 1000;
	//y(j)方向の毛細管作用による浸透距離
	ry = (gap[i][j][X] <= gap[i][j][Z]) ? (gap[i][j][X] / 2) : (gap[i][j][Z] / 2);
	hy = sqrt(((ry * 0.001) * SurfaceTension * cos(ContactAngle) * dt) / (2 * Viscosity)) * 1000;


	/////////////////////////////////////////防染部分(板など)のある方向の探索////////////////////////////////////////////
	int Ncount = 0, Scount = 0, Wcount = 0, Ecount = 0;
	while (p[i][j + Ncount][k] < P_MAX) {
		Ncount++;
		if (j + Ncount > 2 * Nx + 1) break;
	}
	while (p[i][j - Scount][k] < P_MAX) {
		Scount++;
		if (j - Scount < 0) break;
	}
	while (p[i - Wcount][j][k] < P_MAX) {
		Wcount++;
		if (i - Wcount < 0) break;
	}
	while (p[i + Ecount][j][k] < P_MAX) {
		Ecount++;
		if (i + Ecount > 2 * Ny + 1) break;
	}
	double S = 1.5;
	if (Ncount < Scount) {
		DyN *= S;
		DgN *= S;
		DyS /= S;
		DgS /= S;
	}
	else {
		DyN /= S;
		DgN /= S;
		DyS *= S;
		DgS *= S;
	}
	if (Wcount < Ecount) {
		DyW *= S;
		DgW *= S;
		DyE /= S;
		DgE /= S;
	}
	else {
		DyW /= S;
		DgW /= S;
		DyE *= S;
		DgE *= S;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////横方向の糸間の隙間：i = 偶数，j = 奇数//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (i % 2 == 0 && j % 2 == 1) {
		///////////////////////糸への浸透の計算//////////////////////////////
		if (AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
			if (c[i][j][k] > 0) {
				dcN_yarn = dt *  DyN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((yarn[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
				dcS_yarn = dt *  DyS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((yarn[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
				if (dcN_yarn > 0) {
					c[i][j + 1][k] += dcN_yarn;
					dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
					c[i][j][k] -= dcN_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
				if (dcS_yarn > 0) {
					c[i][j - 1][k] += dcS_yarn;
					dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
					c[i][j][k] -= dcS_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
		}
		else if (AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
			//(j + 1)：拡散
			if (c[i][j][k] > 0) {
				dcN_yarn = dt *  DyN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((yarn[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
				if (dcN_yarn > 0) {
					c[i][j + 1][k] += dcN_yarn;
					dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
					c[i][j][k] -= dcN_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			//(j - 1)：浸透, Buras，着目セルに染料溶液があるとき
			if ((dye[i][j][k] + water[i][j][k]) > 0 && (dye[i][j - 1][k] + water[i][j - 1][k]) < capacity[i][j - 1]) {
				elapsedTime = (water[i][j - 1][k] > 0) ? (-(capacity[i][j - 1] / initialVelocity) * log(1 - ((water[i][j - 1][k] + dye[i][j - 1][k]) / capacity[i][j - 1]))) : 0;
				dq = capacity[i][j - 1] * (1 - exp(-(initialVelocity / capacity[i][j - 1]) * dt)) * exp(-(initialVelocity / capacity[i][j - 1]) * elapsedTime);
				if ((dye[i][j][k] + water[i][j][k]) < dq) {
					dq = dye[i][j][k] + water[i][j][k];
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] = dye[i][j][k] = c[i][j][k] = 0;
				}
				else {
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] -= dwater;
					dye[i][j][k] -= ddye;
					c[i][j][k] = dye[i][j][k] / (dye[i][j][k] + water[i][j][k]);
				}
				water[i][j - 1][k] += dwater;
				dye[i][j - 1][k] += ddye;
				c[i][j - 1][k] = dye[i][j - 1][k] / (dye[i][j - 1][k] + water[i][j - 1][k]);
			}
		}
		else if (AdjacentCellStatus[i][j][k][3] == 1 && AdjacentCellStatus[i][j][k][2] == 0) {
			//(j - 1)：拡散
			if (c[i][j][k] > 0) {
				dcS_yarn = dt *  DyS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((yarn[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
				if (dcS_yarn > 0) {
					c[i][j - 1][k] += dcS_yarn;
					dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
					c[i][j][k] -= dcS_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			//(j + 1)：浸透, Buras
			if ((dye[i][j][k] + water[i][j][k]) > 0 && (dye[i][j + 1][k] + water[i][j + 1][k]) < capacity[i][j + 1]) {
				elapsedTime = (water[i][j + 1][k] > 0) ? (-(capacity[i][j + 1] / initialVelocity) * log(1 - ((water[i][j + 1][k] + dye[i][j + 1][k]) / capacity[i][j + 1]))) : 0;
				dq = capacity[i][j + 1] * (1 - exp(-(initialVelocity / capacity[i][j + 1]) * dt)) * exp(-(initialVelocity / capacity[i][j + 1]) * elapsedTime);
				if ((dye[i][j][k] + water[i][j][k]) < dq) {
					dq = dye[i][j][k] + water[i][j][k];
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] = dye[i][j][k] = c[i][j][k] = 0;
				}
				else {
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] -= dwater;
					dye[i][j][k] -= ddye;
					c[i][j][k] = dye[i][j][k] / (dye[i][j][k] + water[i][j][k]);
				}
				water[i][j + 1][k] += dwater;
				dye[i][j + 1][k] += ddye;
				c[i][j + 1][k] = dye[i][j + 1][k] / (dye[i][j + 1][k] + water[i][j + 1][k]);
			}
		}
		else if (AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
			//浸透計算：Burasの式 dq = Q * (1 - exp(-(I / Q) * dt));
			if ((dye[i][j][k] + water[i][j][k]) > 0) {
				elapsedTime = (water[i][j + 1][k] > 0) ? (-(capacity[i][j + 1] / initialVelocity) * log(1 - ((water[i][j + 1][k] + dye[i][j + 1][k]) / capacity[i][j + 1]))) : 0;
				dq = capacity[i][j + 1] * (1 - exp(-(initialVelocity / capacity[i][j + 1]) * dt)) * exp(-(initialVelocity / capacity[i][j + 1]) * elapsedTime);
				dq = ((dye[i][j + 1][k] + water[i][j + 1][k]) < capacity[i][j + 1]) ? dq : 0;
				elapsedTime = (water[i][j - 1][k] > 0) ? (-(capacity[i][j - 1] / initialVelocity) * log(1 - ((water[i][j - 1][k] + dye[i][j - 1][k]) / capacity[i][j - 1]))) : 0;
				dq1 = capacity[i][j - 1] * (1 - exp(-(initialVelocity / capacity[i][j - 1]) * dt)) * exp(-(initialVelocity / capacity[i][j - 1]) * elapsedTime);
				dq1 = ((dye[i][j - 1][k] + water[i][j - 1][k]) < capacity[i][j - 1]) ? dq1 : 0;
				if ((dye[i][j][k] + water[i][j][k]) < (dq + dq1)) {
					dq = dq1 = (dye[i][j][k] + water[i][j][k]) / 2;
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					dwater1 = dq1 * (1 - c[i][j][k]);
					ddye1 = dq1 * c[i][j][k];
					water[i][j][k] = dye[i][j][k] = c[i][j][k] = 0;
				}
				else {
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					dwater1 = dq1 * (1 - c[i][j][k]);
					ddye1 = dq1 * c[i][j][k];
					water[i][j][k] -= (dwater + dwater1);
					dye[i][j][k] -= (ddye + ddye1);
					c[i][j][k] = dye[i][j][k] / (dye[i][j][k] + water[i][j][k]);
				}
				water[i][j + 1][k] += dwater;
				dye[i][j + 1][k] += ddye;
				c[i][j + 1][k] = dye[i][j + 1][k] / (dye[i][j + 1][k] + water[i][j + 1][k]);
				water[i][j - 1][k] += dwater1;
				dye[i][j - 1][k] += ddye1;
				c[i][j - 1][k] = dye[i][j - 1][k] / (dye[i][j - 1][k] + water[i][j - 1][k]);
			}
		}
		///////////////////////隙間への浸透の計算////////////////////////////
		//着目セルの水分量が超過しているか，水分はあるか，全くないかで分類
		if (water[i][j][k] > capacity[i][j]) {
			if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1) {
				//2つの隙間セルは湿潤 → 拡散計算
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0) {
				//乾燥している方(i - 1)へ毛細管作用
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//配列外参照を回避
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				if (hnegativex != 0) {
					ddyewater = dye[i][j][k] + water[i][j][k] - capacity[i][j];
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {//配列外参照を回避
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				//湿潤している方(i + 1)へ拡散計算
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][0] == 0) {
				//乾燥している方(i + 1)へ毛細管作用
				while (hx > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx = 0;
					}
					else hx -= gap[i + hpositivex][j][X];
				}
				if (hpositivex != 0) {
					ddyewater = dye[i][j][k] + water[i][j][k] - capacity[i][j];
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					//浸透するセルへ水分を分配
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				//湿潤している方(i - 1)へ拡散計算
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0){
				//2つの隙間セルは乾燥 → 毛細管作用
				hx1 = hx;
				while (hx1 > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx1 < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx1 = 0;
					}
					else hx1 -= gap[i + hpositivex][j][X];
				}
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 2;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
			}
		}
		else if (water[i][j][k] > 0) {
			if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1) {
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0) {
				//湿潤している方(i + 1)へ拡散計算，乾燥している方(i - 1)へは何もしない
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][0] == 0) {
				//湿潤している方(i - 1)へ拡散計算，乾燥している方(i + 1)へは何もしない
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
					}
				}
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////縦方向の糸間の隙間：i = 奇数，j = 偶数//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (i % 2 == 1 && j % 2 == 0) {
		///////////////////////糸への浸透の計算//////////////////////////////
		if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1) {
			if (c[i][j][k] > 0) {
				dcW_yarn = dt *  DyW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((yarn[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
				dcE_yarn = dt *  DyE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((yarn[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
				if (dcE_yarn > 0) {
					c[i + 1][j][k] += dcE_yarn;
					dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
					c[i][j][k] -= dcE_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

				}
				if (dcW_yarn > 0) {
					c[i - 1][j][k] += dcW_yarn;
					dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
					c[i][j][k] -= dcW_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

				}
			}
		}
		else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0) {
			//(i + 1)：拡散
			if (c[i][j][k] > 0) {
				dcE_yarn = dt *  DyE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((yarn[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
				if (dcE_yarn > 0) {
					c[i + 1][j][k] += dcE_yarn;
					dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
					c[i][j][k] -= dcE_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

				}
			}
			//(i - 1)：浸透, Buras
			if ((dye[i][j][k] + water[i][j][k]) > 0 && (dye[i - 1][j][k] + water[i - 1][j][k]) < capacity[i - 1][j]) {
				elapsedTime = (water[i - 1][j][k] > 0) ? (-(capacity[i - 1][j] / initialVelocity) * log(1 - ((water[i - 1][j][k] + dye[i - 1][j][k]) / capacity[i - 1][j]))) : 0;
				dq = capacity[i - 1][j] * (1 - exp(-(initialVelocity / capacity[i - 1][j]) * dt)) * exp(-(initialVelocity / capacity[i - 1][j]) * elapsedTime);
				if ((dye[i][j][k] + water[i][j][k]) < dq) {
					dq = dye[i][j][k] + water[i][j][k];
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] = dye[i][j][k] = c[i][j][k] = 0;
				}
				else {
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] -= dwater;
					dye[i][j][k] -= ddye;
					c[i][j][k] = dye[i][j][k] / (dye[i][j][k] + water[i][j][k]);
				}
				water[i - 1][j][k] += dwater;
				dye[i - 1][j][k] += ddye;
				c[i - 1][j][k] = dye[i - 1][j][k] / (dye[i - 1][j][k] + water[i - 1][j][k]);
			}
		}
		else if (AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][0] == 0) {
			//(i - 1)：拡散
			if (c[i][j][k] > 0) {
				dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
				if (dcW_yarn > 0) {
					c[i - 1][j][k] += dcW_yarn;
					dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
					c[i][j][k] -= dcW_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

				}
			}
			//(i + 1)：浸透, Buras
			if ((dye[i][j][k] + water[i][j][k]) > 0 && (dye[i + 1][j][k] + water[i + 1][j][k]) < capacity[i + 1][j]) {
				elapsedTime = (water[i + 1][j][k] > 0) ? (-(capacity[i + 1][j] / initialVelocity) * log(1 - ((water[i + 1][j][k] + dye[i + 1][j][k]) / capacity[i + 1][j]))) : 0;
				dq = capacity[i + 1][j] * (1 - exp(-(initialVelocity / capacity[i + 1][j]) * dt)) * exp(-(initialVelocity / capacity[i + 1][j]) * elapsedTime);
				if ((dye[i][j][k] + water[i][j][k]) < dq) {
					dq = dye[i][j][k] + water[i][j][k];
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] = dye[i][j][k] = c[i][j][k] = 0;
				}
				else {
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					water[i][j][k] -= dwater;
					dye[i][j][k] -= ddye;
					c[i][j][k] = dye[i][j][k] / (dye[i][j][k] + water[i][j][k]);
				}
				water[i + 1][j][k] += dwater;
				dye[i + 1][j][k] += ddye;
				c[i + 1][j][k] = dye[i + 1][j][k] / (dye[i + 1][j][k] + water[i + 1][j][k]);
			}
		}
		else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0) {
			if ((dye[i][j][k] + water[i][j][k]) > 0) {
				elapsedTime = (water[i + 1][j][k] > 0) ? (-(capacity[i + 1][j] / initialVelocity) * log(1 - ((water[i + 1][j][k] + dye[i + 1][j][k]) / capacity[i + 1][j]))) : 0;
				dq = capacity[i + 1][j] * (1 - exp(-(initialVelocity / capacity[i + 1][j]) * dt)) * exp(-(initialVelocity / capacity[i + 1][j]) * elapsedTime);
				dq = ((dye[i + 1][j][k] + water[i + 1][j][k]) < capacity[i + 1][j]) ? dq : 0;
				elapsedTime = (water[i - 1][j][k] > 0) ? (-(capacity[i - 1][j] / initialVelocity) * log(1 - ((water[i - 1][j][k] + dye[i - 1][j][k]) / capacity[i - 1][j]))) : 0;
				dq1 = capacity[i - 1][j] * (1 - exp(-(initialVelocity / capacity[i - 1][j]) * dt)) * exp(-(initialVelocity / capacity[i - 1][j]) * elapsedTime);
				dq1 = ((dye[i - 1][j][k] + water[i - 1][j][k]) < capacity[i - 1][j]) ? dq1 : 0;
				if ((dye[i][j][k] + water[i][j][k]) < (dq + dq1)) {
					dq = dq1 = (dye[i][j][k] + water[i][j][k]) / 2;
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					dwater1 = dq1 * (1 - c[i][j][k]);
					ddye1 = dq1 * c[i][j][k];
					water[i][j][k] = dye[i][j][k] = c[i][j][k] = 0;
				}
				else {
					dwater = dq * (1 - c[i][j][k]);
					ddye = dq * c[i][j][k];
					dwater1 = dq1 * (1 - c[i][j][k]);
					ddye1 = dq1 * c[i][j][k];
					water[i][j][k] -= (dwater + dwater1);
					dye[i][j][k] -= (ddye + ddye1);
					c[i][j][k] = dye[i][j][k] / (dye[i][j][k] + water[i][j][k]);
				}
				water[i + 1][j][k] += dwater;
				dye[i + 1][j][k] += ddye;
				c[i + 1][j][k] = dye[i + 1][j][k] / (dye[i + 1][j][k] + water[i + 1][j][k]);
				water[i - 1][j][k] += dwater1;
				dye[i - 1][j][k] += ddye1;
				c[i - 1][j][k] = dye[i - 1][j][k] / (dye[i - 1][j][k] + water[i - 1][j][k]);
			}
		}
		///////////////////////隙間への浸透の計算////////////////////////////
		//着目セルの水分量が超過しているか，水分はあるか，全くないかで分類
		if (water[i][j][k] > capacity[i][j]) {
			if (AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//2つの隙間セルは湿潤 → 拡散計算
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥している方(j - 1)へ毛細管作用
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//配列外参照を回避
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				if (hnegativey != 0) {
					ddyewater = dye[i][j][k] + water[i][j][k] - capacity[i][j];
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤している方(j + 1)へ拡散計算
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][3] == 1 && AdjacentCellStatus[i][j][k][2] == 0) {
				//乾燥している方(j + 1)へ毛細管作用
				while (hy > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy = 0;
					}
					else hy -= gap[i][j + hpositivey][Y];
				}
				if (hpositivey != 0) {
					ddyewater = dye[i][j][k] + water[i][j][k] - capacity[i][j];
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				//湿潤している方(j - 1)へ拡散計算
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				//2つの隙間セルは乾燥 → 毛細管作用
				hy1 = hy;
				while (hy1 > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy1 < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy1 = 0;
					}
					else hy1 -= gap[i][j + hpositivey][Y];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 2;
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
			}
		}
		else if (water[i][j][k] > 0) {
			if (AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤している方(j + 1)へ拡散計算，乾燥している方(j - 1)へは何もしない
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][3] == 1 && AdjacentCellStatus[i][j][k][2] == 0) {
				//湿潤している方(j - 1)へ拡散計算，乾燥している方(j + 1)へは何もしない
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////糸間の隙間（本当に穴）i = 奇数，j = 奇数////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (i % 2 == 1 && j % 2 == 1) {
		if (water[i][j][k] > capacity[i][j]) {
			//近傍の隙間セルが全て湿潤のとき
			if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//拡散計算
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
				//他の層への濃度・水分・染料移動
				ddyewater = dye[i][j][k] + water[i][j][k] - capacity[i][j];
			}

			//近傍の隙間セルのうち，3つが湿潤のとき，湿潤セルへの水分移動はなし
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥(j - 1)へ毛細管作用
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//配列外参照を回避
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 2;
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤(i + 1, i - 1, j + 1)へ拡散			
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//乾燥(j + 1)へ毛細管作用
				while (hy > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy = 0;
					}
					else hy -= gap[i][j + hpositivey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 2;
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				//湿潤(i + 1, i - 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//乾燥(i - 1)へ毛細管作用
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//配列外参照を回避
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 2;
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {//配列外参照を回避
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				//湿潤(i + 1, j + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//乾燥(i + 1)へ毛細管作用
				while (hx > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx = 0;
					}
					else hx -= gap[i + hpositivex][j][X];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 2;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					//浸透するセルへ水分を分配
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				//湿潤(i - 1, j + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//近傍の隙間セルのうち，2つが湿潤のとき，湿潤セルへの水分移動はなし
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥している方(j + 1, j - 1)へ毛細管作用
				hy1 = hy;
				while (hy1 > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy1 < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy1 = 0;
					}
					else hy1 -= gap[i][j + hpositivey][Y];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 3;
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤している方(i + 1, i - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥している方(i - 1, j - 1)へ毛細管作用
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//配列外参照を回避
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//配列外参照を回避
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 3;
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {//配列外参照を回避
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤している方(i + 1, j + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//乾燥している方(i - 1, j + 1)へ毛細管作用
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//配列外参照を回避
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				while (hy > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy = 0;
					}
					else hy -= gap[i][j + hpositivey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 3;
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {//配列外参照を回避
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				//湿潤している方(i + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥している方(i + 1, j - 1)へ毛細管作用
				while (hx > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx = 0;
					}
					else hx -= gap[i + hpositivex][j][X];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//配列外参照を回避
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 3;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					//浸透するセルへ水分を分配
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤している方(i - 1, j + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//乾燥している方(i + 1, j + 1)へ毛細管作用
				while (hx > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx = 0;
					}
					else hx -= gap[i + hpositivex][j][X];
				}
				while (hy > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy = 0;
					}
					else hy -= gap[i][j + hpositivey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 3;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				//湿潤している方(i - 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//乾燥している方(i + 1, i - 1)へ毛細管作用
				hx1 = hx;
				while (hx1 > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx1 < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx1 = 0;
					}
					else hx1 -= gap[i + hpositivex][j][X];
				}
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 3;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				//湿潤している方(j + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//近傍の隙間セルのうち，1つが湿潤のとき，湿潤セルへの水分移動はなし
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥している方(i - 1, j + 1, j - 1)へ毛細管作用
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				hy1 = hy;
				while (hy1 > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy1 < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy1 = 0;
					}
					else hy1 -= gap[i][j + hpositivey][Y];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 4;
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤している方(i + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥している方(i + 1, j + 1, j - 1)へ毛細管作用
				while (hx > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx = 0;
					}
					else hx -= gap[i + hpositivex][j][X];
				}
				hy1 = hy;
				while (hy1 > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy1 < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy1 = 0;
					}
					else hy1 -= gap[i][j + hpositivey][Y];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 4;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤している方(i - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//乾燥している方(i + 1, i - 1, j - 1)へ毛細管作用
				hx1 = hx;
				while (hx1 > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx1 < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx1 = 0;
					}
					else hx1 -= gap[i + hpositivex][j][X];
				}
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 4;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
				//湿潤している方(j + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//乾燥している方(i + 1, i - 1, j + 1)へ毛細管作用
				hx1 = hx;
				while (hx1 > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx1 < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx1 = 0;
					}
					else hx1 -= gap[i + hpositivex][j][X];
				}
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				while (hy > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy = 0;
					}
					else hy -= gap[i][j + hpositivey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 4;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				//湿潤している方(j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//近傍の隙間セルが全て乾燥のとき
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				hx1 = hx;
				while (hx1 > 0) {
					hpositivex++;
					if (i + hpositivex > NN) break;
					if (hx1 < gap[i + hpositivex][j][X]) {
						hpositivex--;
						hx1 = 0;
					}
					else hx1 -= gap[i + hpositivex][j][X];
				}
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				hy1 = hy;
				while (hy1 > 0) {
					hpositivey++;
					if (j + hpositivey > NN) break;
					if (hy1 < gap[i][j + hpositivey][Y]) {
						hpositivey--;
						hy1 = 0;
					}
					else hy1 -= gap[i][j + hpositivey][Y];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;
					if (hy < gap[i][j - hnegativey][Y]) {
						hnegativey--;
						hy = 0;
					}
					else hy -= gap[i][j - hnegativey][Y];
				}
				ddyewater = (dye[i][j][k] + water[i][j][k] - capacity[i][j]) / 5;
				if (hpositivex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				if (hnegativex != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativex * (hnegativex + 1)) * (hnegativex + 1 - l));
						if (i - l >= 0) {
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				if (hpositivey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hpositivey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivey * (hpositivey + 1)) * (hpositivey + 1 - l));
						if (j + l < NN + 2) {
							dye[i][j + l][k] += ddyewater1 * c[i][j][k];
							water[i][j + l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j + l][k] = dye[i][j + l][k] / (dye[i][j + l][k] + water[i][j + l][k]);
						}
					}
				}
				if (hnegativey != 0) {
					dye[i][j][k] -= ddyewater * c[i][j][k];
					water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
					for (int l = 1; l <= hnegativey; l++) {
						ddyewater1 = (ddyewater / (0.5 * hnegativey * (hnegativey + 1)) * (hnegativey + 1 - l));
						if (j - l >= 0) {
							dye[i][j - l][k] += ddyewater1 * c[i][j][k];
							water[i][j - l][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i][j - l][k] = dye[i][j - l][k] / (dye[i][j - l][k] + water[i][j - l][k]);
						}
					}
				}
			}

			//他層の隙間（重なっているところ）へ水分を分配
			dye[i][j][k] -= ddyewater * c[i][j][k];
			water[i][j][k] -= ddyewater * (1 - c[i][j][k]);
			if (k == 0) {
				dye[i][j][k + 1] += ddyewater * c[i][j][k];
				water[i][j][k + 1] += ddyewater * (1 - c[i][j][k]);
				c[i][j][k + 1] = dye[i][j][k + 1] / (dye[i][j][k + 1] + water[i][j][k + 1]);
			}
			else if (k == 1) {
				dye[i][j][k - 1] += ddyewater * c[i][j][k];
				water[i][j][k - 1] += ddyewater * (1 - c[i][j][k]);
				c[i][j][k - 1] = dye[i][j][k - 1] / (dye[i][j][k - 1] + water[i][j][k - 1]);
			}

		}
		else if (water[i][j][k] > 0) {
			//近傍の隙間セルが全て湿潤のとき，水分移動なし，水分分配なし
			if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//近傍の隙間セルのうち，3つが湿潤のとき → 湿潤している方へのみ拡散計算，水分移動なし，水分分配なし
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤(i + 1, i - 1, j + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//湿潤(i + 1, i - 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//湿潤(i + 1, j + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//湿潤(i - 1, j + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//近傍の隙間セルのうち，2つが湿潤のとき → 湿潤している方へのみ拡散計算，水分移動なし，水分分配なし
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤している方(i + 1, i - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤している方(i + 1, j + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//湿潤している方(i + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤している方(i - 1, j + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//湿潤している方(i - 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 1) {
				//湿潤している方(j + 1, j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//近傍の隙間セルのうち，1つが湿潤のとき → 湿潤している方へのみ拡散計算，水分移動なし，水分分配なし
			else if (AdjacentCellStatus[i][j][k][0] == 1 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤している方(i + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//→
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 1 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤している方(i - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//←
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 1 && AdjacentCellStatus[i][j][k][3] == 0) {
				//湿潤している方(j + 1)へ拡散
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//↑
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 1) {
				//湿潤している方(j - 1)へ拡散
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//↓
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//近傍の隙間セルが全て乾燥のとき → なにもない
			else if (AdjacentCellStatus[i][j][k][0] == 0 && AdjacentCellStatus[i][j][k][1] == 0 && AdjacentCellStatus[i][j][k][2] == 0 && AdjacentCellStatus[i][j][k][3] == 0) {}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////糸：i = 偶数，j = 偶数////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else {
		//湿潤している隙間セルへ微量の拡散
		if (c[i][j][k] > 0) {
			DgN = DgS = DgW = DgE = Dg2;
			if (AdjacentCellStatus[i][j][k][0] == 1) {
				dcE = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + yarn[i][j][X]) / 2, 2)));//→
				if (dcE > 0) {
					c[i + 1][j][k] += dcE;
					dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
					c[i][j][k] -= dcE;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (AdjacentCellStatus[i][j][k][1] == 1) {
				dcW = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + yarn[i][j][X]) / 2, 2)));//←
				if (dcW > 0) {
					c[i - 1][j][k] += dcW;
					dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
					c[i][j][k] -= dcW;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (AdjacentCellStatus[i][j][k][2] == 1) {
				dcN = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + yarn[i][j][Y]) / 2, 2)));//↑
				if (dcN > 0) {
					c[i][j + 1][k] += dcN;
					dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
					c[i][j][k] -= dcN;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (AdjacentCellStatus[i][j][k][3] == 1) {
				dcS = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + yarn[i][j][Y]) / 2, 2)));//↓
				if (dcS > 0) {
					c[i][j - 1][k] += dcS;
					dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
					c[i][j][k] -= dcS;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
		}
		//糸の水分・染料溶液量が許容量以上のとき，隙間に糸の水分が染み出す
		dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
		if ((water[i][j][k] + dye[i][j][k]) > capacity[i][j]) {
			//染み出す総量
			ddyewater = (water[i][j][k] + dye[i][j][k] - capacity[i][j]) / 4;
			ddye = ddyewater * c[i][j][k];
			dwater = ddyewater * (1 - c[i][j][k]);

			//近傍セルへの水分分配・濃度計算
			dye[i + 1][j][k] += ddye; water[i + 1][j][k] += dwater;
			c[i + 1][j][k] = dye[i + 1][j][k] / (dye[i + 1][j][k] + water[i + 1][j][k]);

			dye[i - 1][j][k] += ddye; water[i - 1][j][k] += dwater;
			c[i - 1][j][k] = dye[i - 1][j][k] / (dye[i - 1][j][k] + water[i - 1][j][k]);

			dye[i][j + 1][k] += ddye; water[i][j + 1][k] += dwater;
			c[i][j + 1][k] = dye[i][j + 1][k] / (dye[i][j + 1][k] + water[i][j + 1][k]);

			dye[i][j - 1][k] += ddye; water[i][j - 1][k] += dwater;
			c[i][j - 1][k] = dye[i][j - 1][k] / (dye[i][j - 1][k] + water[i][j - 1][k]);

			//糸セルから移動分を差し引く
			water[i][j][k] -= 4 * dwater;
			dye[i][j][k] -= 4 * ddye;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////折り合わされた布への浸透・拡散//////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double dc, d_match = 0.01;
	int match_i = ContactCell_i[i][j][M_OR_V];
	int match_j = ContactCell_j[i][j][M_OR_V];
	if (water[match_i][match_j][k] >= capacity[match_i][match_j] * 1.0) {
		dc = dt * D_match *((c[i][j][k] - c[match_i][match_j][k]) / d_match);
		if (dc > 0) {
			c[match_i][match_j][k] += dc;
			dye[match_i][match_j][k] = (c[match_i][match_j][k] / (1 - c[match_i][match_j][k])) * water[match_i][match_j][k];
			c[i][j][k] -= dc;
			dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////染色終了条件/////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double dThreshold = pow(10, -4) * 8.0;//3.5 * 10^-4が無難，ある一定以上小さくするとなぞの円環が現れる
	if (i % 2 == 0 && j % 2 == 1) {
		if (dcN_yarn > dThreshold || dcS_yarn > dThreshold || dcW_gap > dThreshold || dcE_gap > dThreshold) dCount[i][j][k]++;
	}
	else if (i % 2 == 1 && j % 2 == 0) {
		if (dcN_gap > dThreshold || dcS_gap > dThreshold || dcW_yarn > dThreshold || dcE_yarn > dThreshold) dCount[i][j][k]++;
	}
	else if (i % 2 == 1 && j % 2 == 1) {
		if (dcN_gap > dThreshold || dcS_gap > dThreshold || dcW_gap > dThreshold || dcE_gap > dThreshold) dCount[i][j][k]++;
	}
	double bcThreshold = pow(10, -3);
	//毛細管作用とBurasによる計算が行われたかを計測
	if (dq > bcThreshold) bcCount[i][j][k]++;
	if (dq1 > bcThreshold) bcCount[i][j][k]++;
	if (hpositivex > 0) bcCount[i][j][k]++;
	if (hpositivey > 0) bcCount[i][j][k]++;
	if (hnegativex > 0) bcCount[i][j][k]++;
	if (hnegativey > 0) bcCount[i][j][k]++;

}

//---------------------------------------------------------------------------------------------------
//　　DrawGapPlain
//　　Desc : 描画，平織，隙間の描画なし
//---------------------------------------------------------------------------------------------------
void DrawGapPlain(int i, int j, int k) {

	double gap_coordinate[4][2];//隙間の座標値，第一引数は0:左下，1:右下，2:右上，3:左上，第二引数は0:X，1:Y

	//描画座標値の計算
	//隙間だけわざと大きく設定
	gap_coordinate[0][X] = Gap[i][j][X] - gap[i][j][X];
	gap_coordinate[0][Y] = Gap[i][j][Y] - gap[i][j][Y];
	gap_coordinate[1][X] = gap_coordinate[0][X] + gap[i][j][X] * 3;
	gap_coordinate[1][Y] = gap_coordinate[0][Y];
	gap_coordinate[2][X] = gap_coordinate[1][X];
	gap_coordinate[2][Y] = gap_coordinate[0][Y] + gap[i][j][Y] * 3;
	gap_coordinate[3][X] = gap_coordinate[0][X];
	gap_coordinate[3][Y] = gap_coordinate[2][Y];

	for (int l = 0; l < 4; l++) {
		for (int m = 0; m < 2; m++) {
			gap_coordinate[l][m] *= converting_rate;
			gap_coordinate[l][m] = floor(gap_coordinate[l][m]);
		}
	}

	//布の構造の隙間（縦糸と横糸の格子の隙間）
	double gapcolor = 0;//周囲の糸の平均の色にする
	gapcolor = (dyeDraw[i + 1][j + 1][k] + dyeDraw[i - 1][j - 1][k] + dyeDraw[i + 1][j - 1][k] + dyeDraw[i - 1][j + 1][k]) / 4;
	glColor3d(1 - gapcolor, 1 - gapcolor, 1);
	//glColor3d(1, 0, 0);
	//if (p[i][j][k] == 1.5) glColor3d(0, 0, 0);//木板で挟んだところ
	//if (p[i][j][k] <= 1.0) glColor3d(0, 0, 0);
	if (Folded[i][j] == 2) glColor3d(0, 0, 0);//折れ線
	if (i % 2 == 1 && j % 2 == 1) {
		glPushMatrix();
		glBegin(GL_QUADS);
		glVertex2d(gap_coordinate[0][X], gap_coordinate[0][Y]);//左下
		glVertex2d(gap_coordinate[1][X], gap_coordinate[1][Y]);//右下
		glVertex2d(gap_coordinate[2][X], gap_coordinate[2][Y]);//右上
		glVertex2d(gap_coordinate[3][X], gap_coordinate[3][Y]);//左上
		glEnd();
		glPopMatrix();
	}

}
//---------------------------------------------------------------------------------------------------
//　　DrawYarnPlain
//　　Desc : 描画，平織，隙間の描画なし
//---------------------------------------------------------------------------------------------------
void DrawYarnPlain(int i, int j, int k) {

	double warp_coordinate[4][2];//縦糸の座標値，第一引数は0:左下，1:右下，2:右上，3:左上，第二引数は0:X，1:Y
	double weft_coordinate[4][2];//横糸の座標値，第一引数は0:左下，1:右下，2:右上，3:左上，第二引数は0:X，1:Y	

	//描画座標値の計算
	weft_coordinate[0][X] = weft[i][j][X];
	weft_coordinate[0][Y] = weft[i][j][Y];
	weft_coordinate[1][X] = weft_coordinate[0][X] + weftsize[i][j][X];
	weft_coordinate[1][Y] = weft_coordinate[0][Y];
	weft_coordinate[2][X] = weft_coordinate[1][X];
	weft_coordinate[2][Y] = weft_coordinate[0][Y] + weftsize[i][j][Y];
	weft_coordinate[3][X] = weft_coordinate[0][X];
	weft_coordinate[3][Y] = weft_coordinate[2][Y];

	warp_coordinate[0][X] = warp[i][j][X];
	warp_coordinate[0][Y] = warp[i][j][Y];
	warp_coordinate[1][X] = warp_coordinate[0][X] + warpsize[i][j][X];
	warp_coordinate[1][Y] = warp_coordinate[0][Y];
	warp_coordinate[2][X] = warp_coordinate[1][X];
	warp_coordinate[2][Y] = warp_coordinate[0][Y] + warpsize[i][j][Y];
	warp_coordinate[3][X] = warp_coordinate[0][X];
	warp_coordinate[3][Y] = warp_coordinate[2][Y];

	for (int l = 0; l < 4; l++) {
		for (int m = 0; m < 2; m++) {
			weft_coordinate[l][m] *= converting_rate;
			warp_coordinate[l][m] *= converting_rate;
			weft_coordinate[l][m] = floor(weft_coordinate[l][m]);
			warp_coordinate[l][m] = floor(warp_coordinate[l][m]);
		}
	}

	//糸の描画
	glColor3d(1 - dyeDraw[i][j][k], 1 - dyeDraw[i][j][k], 1);
	//if (p[i][j][k] == 1.5) glColor3d(0, 0, 0);//木板で挟んだところ
	//if (p[i][j][k] <= 1.0) glColor3d(0, 0, 0);
	if (Folded[i][j] == 2) glColor3d(0, 0, 0);//折れ線
	if (k == 0) {//一層目
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//横糸，下			
			glPushMatrix();
			//glColor3d(0, 1, 0);
			glBegin(GL_QUADS);
			glVertex2d(weft_coordinate[0][X], weft_coordinate[0][Y]);//左下
			glVertex2d(weft_coordinate[1][X], weft_coordinate[1][Y]);//右下
			glVertex2d(weft_coordinate[2][X], weft_coordinate[2][Y]);//右上
			glVertex2d(weft_coordinate[3][X], weft_coordinate[3][Y]);//左上
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//縦糸，下
			glPushMatrix();
			//glColor3d(0, 0, 1);
			glBegin(GL_QUADS);
			glVertex2d(warp_coordinate[0][X], warp_coordinate[0][Y]);//左下
			glVertex2d(warp_coordinate[1][X], warp_coordinate[1][Y]);//右下
			glVertex2d(warp_coordinate[2][X], warp_coordinate[2][Y]);//右上
			glVertex2d(warp_coordinate[3][X], warp_coordinate[3][Y]);//左上
			glEnd();
			glPopMatrix();
		}
	}
	else if (k == 1) {//二層目
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//縦糸，上
			glPushMatrix();
			//glColor3d(0, 0, 1);
			glBegin(GL_QUADS);
			glVertex2d(warp_coordinate[0][X], warp_coordinate[0][Y]);//左下
			glVertex2d(warp_coordinate[1][X], warp_coordinate[1][Y]);//右下
			glVertex2d(warp_coordinate[2][X], warp_coordinate[2][Y]);//右上
			glVertex2d(warp_coordinate[3][X], warp_coordinate[3][Y]);//左上
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//横糸，上
			glPushMatrix();
			//glColor3d(0, 1, 0);
			glBegin(GL_QUADS);
			glVertex2d(weft_coordinate[0][X], weft_coordinate[0][Y]);//左下
			glVertex2d(weft_coordinate[1][X], weft_coordinate[1][Y]);//右下
			glVertex2d(weft_coordinate[2][X], weft_coordinate[2][Y]);//右上
			glVertex2d(weft_coordinate[3][X], weft_coordinate[3][Y]);//左上
			glEnd();
			glPopMatrix();
		}
	}

}