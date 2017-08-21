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
//@@PlateBoundaryCal
//@@Desc :	”Â‚Å–hõ‚³‚ê‚½‚Æ‚«‚ÌC”Â‚Ì‹«ŠE•”‚É‚¨‚¯‚éõ—¿‚Ì—¬‚êEZ“§
//---------------------------------------------------------------------------------------------------
void PlateBoundaryCal(int i, int j, int k) {

	//’…–ÚƒZƒ‹‚É‘Î‚µ‚ÄC“™ˆ³—Í‚ÌƒZƒ‹‚ğ’Šo
	int p_equal_i[9] = {}, p_equal_j[9] = {};
	int n = 0;
	for (int l = -1; l <= 1; l++) {
		for (int m = -1; m <= 1; m++) {
			if (p[i][j][k] == p[i + l][j + m][k]) {
				p_equal_i[n] = i + l;
				p_equal_j[n] = j + m;
				n++;
			}
		}
	}
	if (n > 0) {
		double dyesum = 0, csum = 0, watersum = 0;
		for (int l = 0; l < n; l++) {
			dyesum += dye[p_equal_i[l]][p_equal_j[l]][k];
			csum += c[p_equal_i[l]][p_equal_j[l]][k];
			watersum += water[p_equal_i[l]][p_equal_j[l]][k];
		}
		dyesum /= n;
		csum /= n;
		watersum /= n;
		for (int l = 0; l < n; l++) {
			dye[p_equal_i[l]][p_equal_j[l]][k] = dyesum;
			c[p_equal_i[l]][p_equal_j[l]][k] = csum;
			water[p_equal_i[l]][p_equal_j[l]][k] = watersum;
		}
	}

}


//---------------------------------------------------------------------------------------------------
//@@PlainCal
//@@Desc : •½DCõF‚ÌŒvZ‚Ì”»’fC–Ñ×ŠÇ‚È‚Ì‚©ŠgU‚È‚Ì‚©
//---------------------------------------------------------------------------------------------------
void PlainCal(int i, int j, int k) {

	hx = hx1 = hy = hy1 = 0;//’l‚Ì‰Šú‰»
	//ƒZƒ‹‚Ì…•ª—Ê‚ª–O˜a—Ê‚ğ’´‚¦‚½“_‚Åi‚Ç‚ñ‚È‚É­—Ê‚Å‚àj—×‚ÌƒZƒ‹‚É…•ª‚ğó‚¯“n‚·
	hpositivex = hnegativex = hpositivey = hnegativey = 0;
	r = 0;
	dq = dq1 = 0;

	//Œ„ŠÔ ¨ Œ„ŠÔCŒ„ŠÔ ¨ …‚É‚¨‚¯‚éŠgUˆÚ“®—Ê‚ğŒvZ
	dcN = dcS = dcW = dcE = 0;
	dcN_yarn = dcS_yarn = dcW_yarn = dcE_yarn = 0;
	dcN_gap = dcS_gap = dcW_gap = dcE_gap = 0;
	DyN = DyS = DyW = DyE = Dy1;
	DgN = DgS = DgW = DgE = Dg1;
	//x(i)•ûŒü‚Ì–Ñ×ŠÇì—p‚É‚æ‚éZ“§‹——£
	rx = (gap[i][j][Y] <= gap[i][j][Z]) ? (gap[i][j][Y] / 2) : (gap[i][j][Z] / 2);
	hx = sqrt(((rx * 0.001) * SurfaceTension * cos(ContactAngle) * dt) / (2 * Viscosity)) * 1000;
	//y(j)•ûŒü‚Ì–Ñ×ŠÇì—p‚É‚æ‚éZ“§‹——£
	ry = (gap[i][j][X] <= gap[i][j][Z]) ? (gap[i][j][X] / 2) : (gap[i][j][Z] / 2);
	hy = sqrt(((ry * 0.001) * SurfaceTension * cos(ContactAngle) * dt) / (2 * Viscosity)) * 1000;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////‰¡•ûŒü‚Ì…ŠÔ‚ÌŒ„ŠÔFi = ‹ô”Cj = Šï”//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (i % 2 == 0 && j % 2 == 1) {
		///////////////////////…‚Ö‚ÌZ“§‚ÌŒvZ//////////////////////////////
		if (WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
			if (c[i][j][k] > 0) {
				dcN_yarn = dt *  DyN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((yarn[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
				dcS_yarn = dt *  DyS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((yarn[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
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
		else if (WetDryFlag[i][j][k][2] == 1) {
			//(j + 1)FŠgU
			if (c[i][j][k] > 0) {
				dcN_yarn = dt *  DyN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((yarn[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
				if (dcN_yarn > 0) {
					c[i][j + 1][k] += dcN_yarn;
					dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
					c[i][j][k] -= dcN_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			//(j - 1)FZ“§, BurasC’…–ÚƒZƒ‹‚Éõ—¿—n‰t‚ª‚ ‚é‚Æ‚«
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
		else if (WetDryFlag[i][j][k][3] == 1) {
			//(j - 1)FŠgU
			if (c[i][j][k] > 0) {
				dcS_yarn = dt *  DyS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((yarn[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
				if (dcS_yarn > 0) {
					c[i][j - 1][k] += dcS_yarn;
					dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
					c[i][j][k] -= dcS_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			//(j + 1)FZ“§, Buras
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
		else {
			//Z“§ŒvZFBuras‚Ì® dq = Q * (1 - exp(-(I / Q) * dt));
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
		///////////////////////Œ„ŠÔ‚Ö‚ÌZ“§‚ÌŒvZ////////////////////////////
		//’…–ÚƒZƒ‹‚Ì…•ª—Ê‚ª’´‰ß‚µ‚Ä‚¢‚é‚©C…•ª‚Í‚ ‚é‚©C‘S‚­‚È‚¢‚©‚Å•ª—Ş
		if (water[i][j][k] > capacity[i][j]) {
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				//2‚Â‚ÌŒ„ŠÔƒZƒ‹‚Í¼ ¨ ŠgUŒvZ
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i - 1)‚Ö–Ñ×ŠÇì—p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
						if (i - l >= 0) {//”z—ñŠOQÆ‚ğ‰ñ”ğ
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				//¼‚µ‚Ä‚¢‚é•û(i + 1)‚ÖŠgUŒvZ
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i + 1)‚Ö–Ñ×ŠÇì—p
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
					//Z“§‚·‚éƒZƒ‹‚Ö…•ª‚ğ•ª”z
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				//¼‚µ‚Ä‚¢‚é•û(i - 1)‚ÖŠgUŒvZ
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else {
				//2‚Â‚ÌŒ„ŠÔƒZƒ‹‚ÍŠ£‘‡ ¨ –Ñ×ŠÇì—p
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
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i + 1)‚ÖŠgUŒvZCŠ£‘‡‚µ‚Ä‚¢‚é•û(i - 1)‚Ö‚Í‰½‚à‚µ‚È‚¢
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i - 1)‚ÖŠgUŒvZCŠ£‘‡‚µ‚Ä‚¢‚é•û(i + 1)‚Ö‚Í‰½‚à‚µ‚È‚¢
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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
	//////////////////////////////////////////////////c•ûŒü‚Ì…ŠÔ‚ÌŒ„ŠÔFi = Šï”Cj = ‹ô”//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (i % 2 == 1 && j % 2 == 0) {
		///////////////////////…‚Ö‚ÌZ“§‚ÌŒvZ//////////////////////////////
		if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
			if (c[i][j][k] > 0) {
				dcW_yarn = dt *  DyW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((yarn[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
				dcE_yarn = dt *  DyE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((yarn[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
		else if (WetDryFlag[i][j][k][0] == 1) {
			//(i + 1)FŠgU
			if (c[i][j][k] > 0) {
				dcE_yarn = dt *  DyE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((yarn[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
				if (dcE_yarn > 0) {
					c[i + 1][j][k] += dcE_yarn;
					dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
					c[i][j][k] -= dcE_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

				}
			}
			//(i - 1)FZ“§, Buras
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
		else if (WetDryFlag[i][j][k][1] == 1) {
			//(i - 1)FŠgU
			if (c[i][j][k] > 0) {
				dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
				if (dcW_yarn > 0) {
					c[i - 1][j][k] += dcW_yarn;
					dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
					c[i][j][k] -= dcW_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

				}
			}
			//(i + 1)FZ“§, Buras
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
		else {
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
		///////////////////////Œ„ŠÔ‚Ö‚ÌZ“§‚ÌŒvZ////////////////////////////
		//’…–ÚƒZƒ‹‚Ì…•ª—Ê‚ª’´‰ß‚µ‚Ä‚¢‚é‚©C…•ª‚Í‚ ‚é‚©C‘S‚­‚È‚¢‚©‚Å•ª—Ş
		if (water[i][j][k] > capacity[i][j]) {
			if (WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//2‚Â‚ÌŒ„ŠÔƒZƒ‹‚Í¼ ¨ ŠgUŒvZ
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
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
			else if (WetDryFlag[i][j][k][2] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(j - 1)‚Ö–Ñ×ŠÇì—p
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
				//¼‚µ‚Ä‚¢‚é•û(j + 1)‚ÖŠgUŒvZ
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(j + 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(j - 1)‚ÖŠgUŒvZ
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else {
				//2‚Â‚ÌŒ„ŠÔƒZƒ‹‚ÍŠ£‘‡ ¨ –Ñ×ŠÇì—p
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
			if (WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
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
			else if (WetDryFlag[i][j][k][2] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(j + 1)‚ÖŠgUŒvZCŠ£‘‡‚µ‚Ä‚¢‚é•û(j - 1)‚Ö‚Í‰½‚à‚µ‚È‚¢
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][3] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(j - 1)‚ÖŠgUŒvZCŠ£‘‡‚µ‚Ä‚¢‚é•û(j + 1)‚Ö‚Í‰½‚à‚µ‚È‚¢
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
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
	//////////////////////////////////////////////////…ŠÔ‚ÌŒ„ŠÔi–{“–‚ÉŒŠji = Šï”Cj = Šï”////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (i % 2 == 1 && j % 2 == 1) {
		if (water[i][j][k] > capacity[i][j]) {
			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚ª‘S‚Ä¼‚Ì‚Æ‚«
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//ŠgUŒvZ
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
				//‘¼‚Ì‘w‚Ö‚Ì”Z“xE…•ªEõ—¿ˆÚ“®
				ddyewater = dye[i][j][k] + water[i][j][k] - capacity[i][j];
			}

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚Ì‚¤‚¿C3‚Â‚ª¼‚Ì‚Æ‚«C¼ƒZƒ‹‚Ö‚Ì…•ªˆÚ“®‚Í‚È‚µ
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//Š£‘‡(j - 1)‚Ö–Ñ×ŠÇì—p
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
				//¼(i + 1, i - 1, j + 1)‚ÖŠgU			
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡(j + 1)‚Ö–Ñ×ŠÇì—p
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
				//¼(i + 1, i - 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡(i - 1)‚Ö–Ñ×ŠÇì—p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
						if (i - l >= 0) {//”z—ñŠOQÆ‚ğ‰ñ”ğ
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				//¼(i + 1, j + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡(i + 1)‚Ö–Ñ×ŠÇì—p
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
					//Z“§‚·‚éƒZƒ‹‚Ö…•ª‚ğ•ª”z
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				//¼(i - 1, j + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚Ì‚¤‚¿C2‚Â‚ª¼‚Ì‚Æ‚«C¼ƒZƒ‹‚Ö‚Ì…•ªˆÚ“®‚Í‚È‚µ
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(j + 1, j - 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(i + 1, i - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i - 1, j - 1)‚Ö–Ñ×ŠÇì—p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
						if (i - l >= 0) {//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
				//¼‚µ‚Ä‚¢‚é•û(i + 1, j + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i - 1, j + 1)‚Ö–Ñ×ŠÇì—p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
						if (i - l >= 0) {//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
				//¼‚µ‚Ä‚¢‚é•û(i + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i + 1, j - 1)‚Ö–Ñ×ŠÇì—p
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
					if (j - hnegativey <= 0) break;//”z—ñŠOQÆ‚ğ‰ñ”ğ
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
					//Z“§‚·‚éƒZƒ‹‚Ö…•ª‚ğ•ª”z
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
				//¼‚µ‚Ä‚¢‚é•û(i - 1, j + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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
			else if (WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i + 1, j + 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(i - 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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
			else if (WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i + 1, i - 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(j + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
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

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚Ì‚¤‚¿C1‚Â‚ª¼‚Ì‚Æ‚«C¼ƒZƒ‹‚Ö‚Ì…•ªˆÚ“®‚Í‚È‚µ
			else if (WetDryFlag[i][j][k][0] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i - 1, j + 1, j - 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(i + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i + 1, j + 1, j - 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(i - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][2] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i + 1, i - 1, j - 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(j + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][3] == 1) {
				//Š£‘‡‚µ‚Ä‚¢‚é•û(i + 1, i - 1, j + 1)‚Ö–Ñ×ŠÇì—p
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
				//¼‚µ‚Ä‚¢‚é•û(j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚ª‘S‚ÄŠ£‘‡‚Ì‚Æ‚«
			else {
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

			//‘¼‘w‚ÌŒ„ŠÔid‚È‚Á‚Ä‚¢‚é‚Æ‚±‚ëj‚Ö…•ª‚ğ•ª”z
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
			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚ª‘S‚Ä¼‚Ì‚Æ‚«C…•ªˆÚ“®‚È‚µC…•ª•ª”z‚È‚µ
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚Ì‚¤‚¿C3‚Â‚ª¼‚Ì‚Æ‚« ¨ ¼‚µ‚Ä‚¢‚é•û‚Ö‚Ì‚İŠgUŒvZC…•ªˆÚ“®‚È‚µC…•ª•ª”z‚È‚µ
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//¼(i + 1, i - 1, j + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//¼(i + 1, i - 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//¼(i + 1, j + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//¼(i - 1, j + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚Ì‚¤‚¿C2‚Â‚ª¼‚Ì‚Æ‚« ¨ ¼‚µ‚Ä‚¢‚é•û‚Ö‚Ì‚İŠgUŒvZC…•ªˆÚ“®‚È‚µC…•ª•ª”z‚È‚µ
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i + 1, i - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i + 1, j + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
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
			else if (WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i - 1, j + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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
			else if (WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i - 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
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
			else if (WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(j + 1, j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
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

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚Ì‚¤‚¿C1‚Â‚ª¼‚Ì‚Æ‚« ¨ ¼‚µ‚Ä‚¢‚é•û‚Ö‚Ì‚İŠgUŒvZC…•ªˆÚ“®‚È‚µC…•ª•ª”z‚È‚µ
			else if (WetDryFlag[i][j][k][0] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//¨
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(i - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//©
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][2] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(j + 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//ª
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}
			else if (WetDryFlag[i][j][k][3] == 1) {
				//¼‚µ‚Ä‚¢‚é•û(j - 1)‚ÖŠgU
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//«
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];

					}
				}
			}

			//‹ß–T‚ÌŒ„ŠÔƒZƒ‹‚ª‘S‚ÄŠ£‘‡‚Ì‚Æ‚« ¨ ‚È‚É‚à‚È‚¢
			else {}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////…Fi = ‹ô”Cj = ‹ô”////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else {
		//¼‚µ‚Ä‚¢‚éŒ„ŠÔƒZƒ‹‚Ö”÷—Ê‚ÌŠgU
		if (c[i][j][k] > 0) {
			DgN = DgS = DgW = DgE = Dg2;
			if (WetDryFlag[i][j][k][0] == 1) {
				dcE = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + yarn[i][j][X]) / 2, 2)));//¨
				if (dcE > 0) {
					c[i + 1][j][k] += dcE;
					dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
					c[i][j][k] -= dcE;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (WetDryFlag[i][j][k][1] == 1) {
				dcW = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + yarn[i][j][X]) / 2, 2)));//©
				if (dcW > 0) {
					c[i - 1][j][k] += dcW;
					dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
					c[i][j][k] -= dcW;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (WetDryFlag[i][j][k][2] == 1) {
				dcN = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + yarn[i][j][Y]) / 2, 2)));//ª
				if (dcN > 0) {
					c[i][j + 1][k] += dcN;
					dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
					c[i][j][k] -= dcN;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (WetDryFlag[i][j][k][3] == 1) {
				dcS = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + yarn[i][j][Y]) / 2, 2)));//«
				if (dcS > 0) {
					c[i][j - 1][k] += dcS;
					dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
					c[i][j][k] -= dcS;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
		}
		//…‚Ì…•ªEõ—¿—n‰t—Ê‚ª‹–—e—ÊˆÈã‚Ì‚Æ‚«CŒ„ŠÔ‚É…‚Ì…•ª‚ªõ‚İo‚·
		dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
		if ((water[i][j][k] + dye[i][j][k]) > capacity[i][j]) {
			//õ‚İo‚·‘—Ê
			ddyewater = (water[i][j][k] + dye[i][j][k] - capacity[i][j]) / 4;
			ddye = ddyewater * c[i][j][k];
			dwater = ddyewater * (1 - c[i][j][k]);

			//‹ß–TƒZƒ‹‚Ö‚Ì…•ª•ª”zE”Z“xŒvZ
			dye[i + 1][j][k] += ddye; water[i + 1][j][k] += dwater;
			c[i + 1][j][k] = dye[i + 1][j][k] / (dye[i + 1][j][k] + water[i + 1][j][k]);

			dye[i - 1][j][k] += ddye; water[i - 1][j][k] += dwater;
			c[i - 1][j][k] = dye[i - 1][j][k] / (dye[i - 1][j][k] + water[i - 1][j][k]);

			dye[i][j + 1][k] += ddye; water[i][j + 1][k] += dwater;
			c[i][j + 1][k] = dye[i][j + 1][k] / (dye[i][j + 1][k] + water[i][j + 1][k]);

			dye[i][j - 1][k] += ddye; water[i][j - 1][k] += dwater;
			c[i][j - 1][k] = dye[i][j - 1][k] / (dye[i][j - 1][k] + water[i][j - 1][k]);

			//…ƒZƒ‹‚©‚çˆÚ“®•ª‚ğ·‚µˆø‚­
			water[i][j][k] -= 4 * dwater;
			dye[i][j][k] -= 4 * ddye;
		}
	}


	//õFI—¹ğŒ
	double dThreshold = pow(10, -4) * 3.5;//10^-4‚ª–³“ïC‚ ‚éˆê’èˆÈã¬‚³‚­‚·‚é‚Æ‚È‚¼‚Ì‰~ŠÂ‚ªŒ»‚ê‚é
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
	//–Ñ×ŠÇì—p‚ÆBuras‚É‚æ‚éŒvZ‚ªs‚í‚ê‚½‚©‚ğŒv‘ª
	if (dq > bcThreshold) bcCount[i][j][k]++;
	if (dq1 > bcThreshold) bcCount[i][j][k]++;
	if (hpositivex > 0) bcCount[i][j][k]++;
	if (hpositivey > 0) bcCount[i][j][k]++;
	if (hnegativex > 0) bcCount[i][j][k]++;
	if (hnegativey > 0) bcCount[i][j][k]++;

}

//---------------------------------------------------------------------------------------------------
//@@DrawGapPlain
//@@Desc : •`‰æC•½DCŒ„ŠÔ‚Ì•`‰æ‚È‚µ
//---------------------------------------------------------------------------------------------------
void DrawGapPlain(int i, int j, int k) {

	double gap_coordinate[4][2];//Œ„ŠÔ‚ÌÀ•W’lC‘æˆêˆø”‚Í0:¶‰ºC1:‰E‰ºC2:‰EãC3:¶ãC‘æ“ñˆø”‚Í0:XC1:Y

	//•`‰æÀ•W’l‚ÌŒvZ
	//Œ„ŠÔ‚¾‚¯‚í‚´‚Æ‘å‚«‚­İ’è
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

	//•z‚Ì\‘¢‚ÌŒ„ŠÔic…‚Æ‰¡…‚ÌŠiq‚ÌŒ„ŠÔj
	double gapcolor = 0;//üˆÍ‚Ì…‚Ì•½‹Ï‚ÌF‚É‚·‚é
	gapcolor = (dyeDraw[i + 1][j + 1][k] + dyeDraw[i - 1][j - 1][k] + dyeDraw[i + 1][j - 1][k] + dyeDraw[i - 1][j + 1][k]) / 4;
	glColor3d(1 - gapcolor, 1 - gapcolor, 1);
	//glColor3d(1, 0, 0);
	if (p[i][j][k] == 1.5) glColor3d(0.6, 0.2, 0.2);//–Ø”Â‚Å‹²‚ñ‚¾‚Æ‚±‚ë
	//if (p[i][j][k] <= 1.0) glColor3d(0, 0, 0);
	if (i % 2 == 1 && j % 2 == 1) {
		glPushMatrix();
		glBegin(GL_QUADS);
		glVertex2d(gap_coordinate[0][X], gap_coordinate[0][Y]);//¶‰º
		glVertex2d(gap_coordinate[1][X], gap_coordinate[1][Y]);//‰E‰º
		glVertex2d(gap_coordinate[2][X], gap_coordinate[2][Y]);//‰Eã
		glVertex2d(gap_coordinate[3][X], gap_coordinate[3][Y]);//¶ã
		glEnd();
		glPopMatrix();
	}

}
//---------------------------------------------------------------------------------------------------
//@@DrawYarnPlain
//@@Desc : •`‰æC•½DCŒ„ŠÔ‚Ì•`‰æ‚È‚µ
//---------------------------------------------------------------------------------------------------
void DrawYarnPlain(int i, int j, int k) {

	double warp_coordinate[4][2];//c…‚ÌÀ•W’lC‘æˆêˆø”‚Í0:¶‰ºC1:‰E‰ºC2:‰EãC3:¶ãC‘æ“ñˆø”‚Í0:XC1:Y
	double weft_coordinate[4][2];//‰¡…‚ÌÀ•W’lC‘æˆêˆø”‚Í0:¶‰ºC1:‰E‰ºC2:‰EãC3:¶ãC‘æ“ñˆø”‚Í0:XC1:Y	

	//•`‰æÀ•W’l‚ÌŒvZ
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

	//…‚Ì•`‰æ
	glColor3d(1 - dyeDraw[i][j][k], 1 - dyeDraw[i][j][k], 1);
	if (p[i][j][k] == 1.5) glColor3d(0.6, 0.2, 0.2);//–Ø”Â‚Å‹²‚ñ‚¾‚Æ‚±‚ë
	//if (p[i][j][k] <= 1.0) glColor3d(0, 0, 0);
	if (k == 0) {//ˆê‘w–Ú
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//‰¡…C‰º			
			glPushMatrix();
			//glColor3d(0, 1, 0);
			glBegin(GL_QUADS);
			glVertex2d(weft_coordinate[0][X], weft_coordinate[0][Y]);//¶‰º
			glVertex2d(weft_coordinate[1][X], weft_coordinate[1][Y]);//‰E‰º
			glVertex2d(weft_coordinate[2][X], weft_coordinate[2][Y]);//‰Eã
			glVertex2d(weft_coordinate[3][X], weft_coordinate[3][Y]);//¶ã
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//c…C‰º
			glPushMatrix();
			//glColor3d(0, 0, 1);
			glBegin(GL_QUADS);
			glVertex2d(warp_coordinate[0][X], warp_coordinate[0][Y]);//¶‰º
			glVertex2d(warp_coordinate[1][X], warp_coordinate[1][Y]);//‰E‰º
			glVertex2d(warp_coordinate[2][X], warp_coordinate[2][Y]);//‰Eã
			glVertex2d(warp_coordinate[3][X], warp_coordinate[3][Y]);//¶ã
			glEnd();
			glPopMatrix();
		}
	}
	else if (k == 1) {//“ñ‘w–Ú
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//c…Cã
			glPushMatrix();
			//glColor3d(0, 0, 1);
			glBegin(GL_QUADS);
			glVertex2d(warp_coordinate[0][X], warp_coordinate[0][Y]);//¶‰º
			glVertex2d(warp_coordinate[1][X], warp_coordinate[1][Y]);//‰E‰º
			glVertex2d(warp_coordinate[2][X], warp_coordinate[2][Y]);//‰Eã
			glVertex2d(warp_coordinate[3][X], warp_coordinate[3][Y]);//¶ã
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//‰¡…Cã
			glPushMatrix();
			//glColor3d(0, 1, 0);
			glBegin(GL_QUADS);
			glVertex2d(weft_coordinate[0][X], weft_coordinate[0][Y]);//¶‰º
			glVertex2d(weft_coordinate[1][X], weft_coordinate[1][Y]);//‰E‰º
			glVertex2d(weft_coordinate[2][X], weft_coordinate[2][Y]);//‰Eã
			glVertex2d(weft_coordinate[3][X], weft_coordinate[3][Y]);//¶ã
			glEnd();
			glPopMatrix();
		}
	}

}