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
//�@�@PlainCal
//�@�@Desc : ���D�C���F�̌v�Z�̔��f�C�э׊ǂȂ̂��g�U�Ȃ̂�
//---------------------------------------------------------------------------------------------------
void PlainCal(int i, int j, int k) {
	
	hx = hx1 = hy = hy1 = 0;//�l�̏�����
	//�Z���̐����ʂ��O�a�ʂ𒴂������_�Łi�ǂ�Ȃɏ��ʂł��j�ׂ̃Z���ɐ������󂯓n��
	hpositivex = hnegativex = hpositivey = hnegativey = 0;
	r = 0;
	dq = dq1 = 0;
	//for (int x = 0; x < 5; x++) WetDryFlag[i][j][k][x] = 0;//M2��5��ł����킹�ł͂����L���ɂ����Ă��܂��Ă������߁C���ꂢ�Ȍ��ʂ��o��

	//���� �� ���ԁC���� �� ���ɂ�����g�U�ړ��ʂ��v�Z
	dcN = dcS = dcW = dcE = 0;
	dcN_yarn = dcS_yarn = dcW_yarn = dcE_yarn = 0;
	dcN_gap = dcS_gap = dcW_gap = dcE_gap = 0;
	DyN = DyS = DyW = DyE = Dy1;
	DgN = DgS = DgW = DgE = Dg1;
	//x(i)�����̖э׊Ǎ�p�ɂ��Z������
	rx = (gap[i][j][Y] <= gap[i][j][Z]) ? (gap[i][j][Y] / 2) : (gap[i][j][Z] / 2);
	hx = sqrt(((rx * 0.001) * SurfaceTension * cos(ContactAngle) * dt) / (2 * Viscosity)) * 1000;
	//y(j)�����̖э׊Ǎ�p�ɂ��Z������
	ry = (gap[i][j][X] <= gap[i][j][Z]) ? (gap[i][j][X] / 2) : (gap[i][j][Z] / 2);
	hy = sqrt(((ry * 0.001) * SurfaceTension * cos(ContactAngle) * dt) / (2 * Viscosity)) * 1000;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////�������̎��Ԃ̌��ԁFi = �����Cj = �//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (i % 2 == 0 && j % 2 == 1) {
		///////////////////////���ւ̐Z���̌v�Z//////////////////////////////
		if (WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
			if (c[i][j][k] > 0) {
				dcN_yarn = dt *  DyN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((yarn[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
				dcS_yarn = dt *  DyS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((yarn[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
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
			//(j + 1)�F�g�U
			if (c[i][j][k] > 0) {
				dcN_yarn = dt *  DyN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((yarn[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
				if (dcN_yarn > 0) {
					c[i][j + 1][k] += dcN_yarn;
					dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
					c[i][j][k] -= dcN_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			//(j - 1)�F�Z��, Buras�C���ڃZ���ɐ����n�t������Ƃ�
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
			//(j - 1)�F�g�U
			if (c[i][j][k] > 0) {
				dcS_yarn = dt *  DyS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((yarn[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
				if (dcS_yarn > 0) {
					c[i][j - 1][k] += dcS_yarn;
					dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
					c[i][j][k] -= dcS_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			//(j + 1)�F�Z��, Buras
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
			//�Z���v�Z�FBuras�̎� dq = Q * (1 - exp(-(I / Q) * dt));
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
		///////////////////////���Ԃւ̐Z���̌v�Z////////////////////////////
		//���ڃZ���̐����ʂ����߂��Ă��邩�C�����͂��邩�C�S���Ȃ����ŕ���
		if (water[i][j][k] > capacity[i][j]) {
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				//2�̌��ԃZ���͎��� �� �g�U�v�Z
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i - 1)�֖э׊Ǎ�p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//�z��O�Q�Ƃ����
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
						if (i - l >= 0) {//�z��O�Q�Ƃ����
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				//�������Ă����(i + 1)�֊g�U�v�Z
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//�������Ă����(i + 1)�֖э׊Ǎ�p
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
					//�Z������Z���֐����𕪔z
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				//�������Ă����(i - 1)�֊g�U�v�Z
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else {
					//2�̌��ԃZ���͊��� �� �э׊Ǎ�p
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
		else if (water[i][j][k] > 0){
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i + 1)�֊g�U�v�Z�C�������Ă����(i - 1)�ւ͉������Ȃ�
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//�������Ă����(i - 1)�֊g�U�v�Z�C�������Ă����(i + 1)�ւ͉������Ȃ�
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
	//////////////////////////////////////////////////�c�����̎��Ԃ̌��ԁFi = ��Cj = ����//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (i % 2 == 1 && j % 2 == 0) {
		///////////////////////���ւ̐Z���̌v�Z//////////////////////////////
		if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
			if (c[i][j][k] > 0) {
				dcW_yarn = dt *  DyW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((yarn[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
				dcE_yarn = dt *  DyE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((yarn[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
			//(i + 1)�F�g�U
			if (c[i][j][k] > 0) {
				dcE_yarn = dt *  DyE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((yarn[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
				if (dcE_yarn > 0) {
					c[i + 1][j][k] += dcE_yarn;
					dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
					c[i][j][k] -= dcE_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
					
				}
			}
			//(i - 1)�F�Z��, Buras
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
			//(i - 1)�F�g�U
			if (c[i][j][k] > 0) {
				dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
				if (dcW_yarn > 0) {
					c[i - 1][j][k] += dcW_yarn;
					dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
					c[i][j][k] -= dcW_yarn;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
					
				}
			}
			//(i + 1)�F�Z��, Buras
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
		else{
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
		///////////////////////���Ԃւ̐Z���̌v�Z////////////////////////////
		//���ڃZ���̐����ʂ����߂��Ă��邩�C�����͂��邩�C�S���Ȃ����ŕ���
		if (water[i][j][k] > capacity[i][j]) {
			if (WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//2�̌��ԃZ���͎��� �� �g�U�v�Z
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
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
					//�������Ă����(j - 1)�֖э׊Ǎ�p
					while (hy > 0) {
						hnegativey++;
						if (j - hnegativey <= 0) break;//�z��O�Q�Ƃ����
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
					//�������Ă����(j + 1)�֊g�U�v�Z
					if (c[i][j][k] > 0) {
						dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
						if (dcN_gap > 0) {
							c[i][j + 1][k] += dcN_gap;
							dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
							c[i][j][k] -= dcN_gap;
							dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
							
						}
					}
				}
			else if (WetDryFlag[i][j][k][3] == 1) {
					//�������Ă����(j + 1)�֖э׊Ǎ�p
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
					//�������Ă����(j - 1)�֊g�U�v�Z
					if (c[i][j][k] > 0) {
						dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
						if (dcS_gap > 0) {
							c[i][j - 1][k] += dcS_gap;
							dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
							c[i][j][k] -= dcS_gap;
							dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
							
						}
					}
				}
			else {
					//2�̌��ԃZ���͊��� �� �э׊Ǎ�p
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
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
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
				//�������Ă����(j + 1)�֊g�U�v�Z�C�������Ă����(j - 1)�ւ͉������Ȃ�
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][3] == 1) {
				//�������Ă����(j - 1)�֊g�U�v�Z�C�������Ă����(j + 1)�ւ͉������Ȃ�
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
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
	//////////////////////////////////////////////////���Ԃ̌��ԁi�{���Ɍ��ji = ��Cj = �////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (i % 2 == 1 && j % 2 == 1) {
		if (water[i][j][k] > capacity[i][j]) {
			//�ߖT�̌��ԃZ�����S�Ď����̂Ƃ�
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				//�g�U�v�Z
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//���̑w�ւ̔Z�x�E�����E�����ړ�
				ddyewater = dye[i][j][k] + water[i][j][k] - capacity[i][j];
			}
			
			//�ߖT�̌��ԃZ���̂����C3�������̂Ƃ��C�����Z���ւ̐����ړ��͂Ȃ�
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//����(j - 1)�֖э׊Ǎ�p
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//�z��O�Q�Ƃ����
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
				//����(i + 1, i - 1, j + 1)�֊g�U			
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//����(j + 1)�֖э׊Ǎ�p
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
				//����(i + 1, i - 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//����(i - 1)�֖э׊Ǎ�p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//�z��O�Q�Ƃ����
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
						if (i - l >= 0) {//�z��O�Q�Ƃ����
							dye[i - l][j][k] += ddyewater1 * c[i][j][k];
							water[i - l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i - l][j][k] = dye[i - l][j][k] / (dye[i - l][j][k] + water[i - l][j][k]);
						}
					}
				}
				//����(i + 1, j + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//����(i + 1)�֖э׊Ǎ�p
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
					//�Z������Z���֐����𕪔z
					for (int l = 1; l <= hpositivex; l++) {
						ddyewater1 = (ddyewater / (0.5 * hpositivex * (hpositivex + 1)) * (hpositivex + 1 - l));
						if (i + l < NN + 2) {
							dye[i + l][j][k] += ddyewater1 * c[i][j][k];
							water[i + l][j][k] += ddyewater1 * (1 - c[i][j][k]);
							c[i + l][j][k] = dye[i + l][j][k] / (dye[i + l][j][k] + water[i + l][j][k]);
						}
					}
				}
				//����(i - 1, j + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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

			//�ߖT�̌��ԃZ���̂����C2�������̂Ƃ��C�����Z���ւ̐����ړ��͂Ȃ�
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				//�������Ă����(j + 1, j - 1)�֖э׊Ǎ�p
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
				//�������Ă����(i + 1, i - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i - 1, j - 1)�֖э׊Ǎ�p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//�z��O�Q�Ƃ����
					if (hx < gap[i - hnegativex][j][X]) {
						hnegativex--;
						hx = 0;
					}
					else hx -= gap[i - hnegativex][j][X];
				}
				while (hy > 0) {
					hnegativey++;
					if (j - hnegativey <= 0) break;//�z��O�Q�Ƃ����
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
						if (i - l >= 0) {//�z��O�Q�Ƃ����
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
				//�������Ă����(i + 1, j + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i - 1, j + 1)�֖э׊Ǎ�p
				while (hx > 0) {
					hnegativex++;
					if (i - hnegativex <= 0) break;//�z��O�Q�Ƃ����
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
						if (i - l >= 0) {//�z��O�Q�Ƃ����
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
				//�������Ă����(i + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i + 1, j - 1)�֖э׊Ǎ�p
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
					if (j - hnegativey <= 0) break;//�z��O�Q�Ƃ����
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
					//�Z������Z���֐����𕪔z
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
				//�������Ă����(i - 1, j + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i + 1, j + 1)�֖э׊Ǎ�p
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
				//�������Ă����(i - 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i + 1, i - 1)�֖э׊Ǎ�p
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
				//�������Ă����(j + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
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

			//�ߖT�̌��ԃZ���̂����C1�������̂Ƃ��C�����Z���ւ̐����ړ��͂Ȃ�
			else if (WetDryFlag[i][j][k][0] == 1) {
				//�������Ă����(i - 1, j + 1, j - 1)�֖э׊Ǎ�p
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
				//�������Ă����(i + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//�������Ă����(i + 1, j + 1, j - 1)�֖э׊Ǎ�p
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
				//�������Ă����(i - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][2] == 1) {
				//�������Ă����(i + 1, i - 1, j - 1)�֖э׊Ǎ�p
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
				//�������Ă����(j + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][3] == 1) {
				//�������Ă����(i + 1, i - 1, j + 1)�֖э׊Ǎ�p
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
				//�������Ă����(j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}

			//�ߖT�̌��ԃZ�����S�Ċ����̂Ƃ�
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

			//���w�̌��ԁi�d�Ȃ��Ă���Ƃ���j�֐����𕪔z
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
			//�ߖT�̌��ԃZ�����S�Ď����̂Ƃ��C�����ړ��Ȃ��C�������z�Ȃ�
			if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1 && WetDryFlag[i][j][k][3] == 1) {
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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

			//�ߖT�̌��ԃZ���̂����C3�������̂Ƃ� �� �������Ă�����ւ̂݊g�U�v�Z�C�����ړ��Ȃ��C�������z�Ȃ�
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1 && WetDryFlag[i][j][k][2] == 1) {
				//����(i + 1, i - 1, j + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//����(i + 1, i - 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//����(i + 1, j + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//����(i - 1, j + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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

			//�ߖT�̌��ԃZ���̂����C2�������̂Ƃ� �� �������Ă�����ւ̂݊g�U�v�Z�C�����ړ��Ȃ��C�������z�Ȃ�
			else if (WetDryFlag[i][j][k][0] == 1 && WetDryFlag[i][j][k][1] == 1) {
				//�������Ă����(i + 1, i - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i + 1, j + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i - 1, j + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(i - 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
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
				//�������Ă����(j + 1, j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
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

			//�ߖT�̌��ԃZ���̂����C1�������̂Ƃ� �� �������Ă�����ւ̂݊g�U�v�Z�C�����ړ��Ȃ��C�������z�Ȃ�
			else if (WetDryFlag[i][j][k][0] == 1) {
				//�������Ă����(i + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcE_gap = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					if (dcE_gap > 0) {
						c[i + 1][j][k] += dcE_gap;
						dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
						c[i][j][k] -= dcE_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][1] == 1) {
				//�������Ă����(i - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcW_gap = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + gap[i][j][X]) / 2, 2)));//��
					if (dcW_gap > 0) {
						c[i - 1][j][k] += dcW_gap;
						dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
						c[i][j][k] -= dcW_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][2] == 1) {
				//�������Ă����(j + 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcN_gap = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					if (dcN_gap > 0) {
						c[i][j + 1][k] += dcN_gap;
						dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
						c[i][j][k] -= dcN_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}
			else if (WetDryFlag[i][j][k][3] == 1) {
				//�������Ă����(j - 1)�֊g�U
				if (c[i][j][k] > 0) {
					dcS_gap = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + gap[i][j][Y]) / 2, 2)));//��
					if (dcS_gap > 0) {
						c[i][j - 1][k] += dcS_gap;
						dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
						c[i][j][k] -= dcS_gap;
						dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
						
					}
				}
			}

			//�ߖT�̌��ԃZ�����S�Ċ����̂Ƃ� �� �Ȃɂ��Ȃ�
			else {}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////���Fi = �����Cj = ����////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else {
		//�������Ă��錄�ԃZ���֔��ʂ̊g�U
		if (c[i][j][k] > 0) {
			DgN = DgS = DgW = DgE = Dg2;
			if (WetDryFlag[i][j][k][0] == 1) {
				dcE = dt *  DgE * ((c[i][j][k] - c[i + 1][j][k]) / (pow((gap[i + 1][j][X] + yarn[i][j][X]) / 2, 2)));//��
				if (dcE > 0) {
					c[i + 1][j][k] += dcE;
					dye[i + 1][j][k] = (c[i + 1][j][k] / (1 - c[i + 1][j][k])) * water[i + 1][j][k];
					c[i][j][k] -= dcE;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (WetDryFlag[i][j][k][1] == 1) {
				dcW = dt *  DgW * ((c[i][j][k] - c[i - 1][j][k]) / (pow((gap[i - 1][j][X] + yarn[i][j][X]) / 2, 2)));//��
				if (dcW > 0) {
					c[i - 1][j][k] += dcW;
					dye[i - 1][j][k] = (c[i - 1][j][k] / (1 - c[i - 1][j][k])) * water[i - 1][j][k];
					c[i][j][k] -= dcW;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (WetDryFlag[i][j][k][2] == 1) {
				dcN = dt *  DgN * ((c[i][j][k] - c[i][j + 1][k]) / (pow((gap[i][j + 1][Y] + yarn[i][j][Y]) / 2, 2)));//��
				if (dcN > 0) {
					c[i][j + 1][k] += dcN;
					dye[i][j + 1][k] = (c[i][j + 1][k] / (1 - c[i][j + 1][k])) * water[i][j + 1][k];
					c[i][j][k] -= dcN;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
			if (WetDryFlag[i][j][k][3] == 1) {
				dcS = dt *  DgS * ((c[i][j][k] - c[i][j - 1][k]) / (pow((gap[i][j - 1][Y] + yarn[i][j][Y]) / 2, 2)));//��
				if (dcS > 0) {
					c[i][j - 1][k] += dcS;
					dye[i][j - 1][k] = (c[i][j - 1][k] / (1 - c[i][j - 1][k])) * water[i][j - 1][k];
					c[i][j][k] -= dcS;
					dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
				}
			}
		}
		//���̐����E�����n�t�ʂ����e�ʈȏ�̂Ƃ��C���ԂɎ��̐��������ݏo��
		dye[i][j][k] = (c[i][j][k] / (1 - c[i][j][k])) * water[i][j][k];
		if ((water[i][j][k] + dye[i][j][k]) > capacity[i][j]) {
			//���ݏo������
			ddyewater = (water[i][j][k] + dye[i][j][k] - capacity[i][j]) / 4;
			ddye = ddyewater * c[i][j][k];
			dwater = ddyewater * (1 - c[i][j][k]);

			//�ߖT�Z���ւ̐������z�E�Z�x�v�Z
			dye[i + 1][j][k] += ddye; water[i + 1][j][k] += dwater;
			c[i + 1][j][k] = dye[i + 1][j][k] / (dye[i + 1][j][k] + water[i + 1][j][k]);

			dye[i - 1][j][k] += ddye; water[i - 1][j][k] += dwater;
			c[i - 1][j][k] = dye[i - 1][j][k] / (dye[i - 1][j][k] + water[i - 1][j][k]);

			dye[i][j + 1][k] += ddye; water[i][j + 1][k] += dwater;
			c[i][j + 1][k] = dye[i][j + 1][k] / (dye[i][j + 1][k] + water[i][j + 1][k]);

			dye[i][j - 1][k] += ddye; water[i][j - 1][k] += dwater;
			c[i][j - 1][k] = dye[i][j - 1][k] / (dye[i][j - 1][k] + water[i][j - 1][k]);

			//���Z������ړ�������������
			water[i][j][k] -= 4 * dwater;
			dye[i][j][k] -= 4 * ddye;
		}
	}


	//���F�I������
	double dThreshold = pow(10, -4) * 3.5;//10^-4������C������ȏ㏬��������ƂȂ��̉~�������
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
	//�э׊Ǎ�p��Buras�ɂ��v�Z���s��ꂽ�����v��
	if (dq > bcThreshold) bcCount[i][j][k]++;
	if (dq1 > bcThreshold) bcCount[i][j][k]++;
	if (hpositivex > 0) bcCount[i][j][k]++;
	if (hpositivey > 0) bcCount[i][j][k]++;
	if (hnegativex > 0) bcCount[i][j][k]++;
	if (hnegativey > 0) bcCount[i][j][k]++;

}

//---------------------------------------------------------------------------------------------------
//�@�@DrawPlain
//�@�@Desc : �`��C���D�C���Ԃ̕`��Ȃ�
//---------------------------------------------------------------------------------------------------
void DrawPlain(int i, int j, int k) {

	I = i / 2; J = j / 2;
	double X1, Y1, X2, Y2, X3, Y3, X4, Y4;
	//�z�̍\���̌��ԁi�c���Ɖ����̊i�q�̌��ԁj
	double gapcolor = 0;
	if (i % 2 == 1 && j % 2 == 1) {
		//���͂̎��̕��ς̐F�ɂ���C4
		gapcolor = (dyeDraw[i + 1][j + 1][k] + dyeDraw[i - 1][j - 1][k] + dyeDraw[i + 1][j - 1][k] + dyeDraw[i - 1][j + 1][k]) / 4;
		//gapcolor = (c[i + 1][j + 1][k] + c[i - 1][j - 1][k] + c[i + 1][j - 1][k] + c[i - 1][j + 1][k]) / 4;
		//gapcolor = 100 * (water[i + 1][j + 1][k] + water[i - 1][j - 1][k] + water[i + 1][j - 1][k] + water[i - 1][j + 1][k]) / 4;
		
		X1 = converting_rate * Gap[i][j][X]; Y1 = converting_rate * Gap[i][j][Y];
		X2 = converting_rate * (Gap[i][j][X] + gap[i][j][X]); Y2 = converting_rate * Gap[i][j][Y];
		X3 = converting_rate * (Gap[i][j][X] + gap[i][j][X]); Y3 = converting_rate * (Gap[i][j][Y] + gap[i][j][Y]);
		X4 = converting_rate * Gap[i][j][X]; Y4 = converting_rate * (Gap[i][j][Y] + gap[i][j][Y]);
		/*X1 = floor(X1); Y1 = floor(Y1);
		X2 = floor(X2); Y2 = floor(Y2);
		X3 = floor(X3); Y3 = floor(Y3);
		X4 = floor(X4); Y4 = floor(Y4);*/

		glPushMatrix();
		glColor3d(1 - gapcolor, 1 - gapcolor, 1);
		glBegin(GL_QUADS);//(I,J)���W�ɂ����āC1�̎��Z���ɂ��C���̉E��Ɍ��Ԃ�`��
		glVertex2d(X1, Y1);//�����`�̍����̊p
		glVertex2d(X2, Y2);//�E��
		glVertex2d(X3, Y3);//�E��
		glVertex2d(X4, Y4);//����
		
		glEnd();
		glPopMatrix();
	}

	//�`��̏��Ԃŏc���Ɖ����̏d�Ȃ��\�����邵���Ȃ�
	//glVertex3d��z�̒l��ω������č������Ă����܂������Ȃ�
	if (k == 0) {//��w��
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//�����C��			
			X1 = converting_rate * weft[i][j][X]; Y1 = converting_rate * weft[i][j][Y];
			X2 = converting_rate * (weft[i][j][X] + 0.5 * gap[i - 1][j][X] + yarn[i][j][X] + 0.5 * gap[i + 1][j][X]); Y2 = converting_rate * weft[i][j][Y];
			X3 = converting_rate * (weft[i][j][X] + 0.5 * gap[i - 1][j][X] + yarn[i][j][X] + 0.5 * gap[i + 1][j][X]); Y3 = converting_rate * (weft[i][j][Y] + yarn[i][j][Y]);
			X4 = converting_rate * weft[i][j][X]; Y4 = converting_rate * (weft[i][j][Y] + yarn[i][j][Y]);
			/*X1 = floor(X1); Y1 = floor(Y1);
			X2 = floor(X2); Y2 = floor(Y2);
			X3 = floor(X3); Y3 = floor(Y3);
			X4 = floor(X4); Y4 = floor(Y4);*/
			glPushMatrix();
			glColor3d(1 - dyeDraw[i][j][k], 1 - dyeDraw[i][j][k], 1);
			//glColor3d(1 - c[i][j][k], 1 - c[i][j][k], 1);
			//glColor3d(1 - 100 * water[i][j][k], 1 - 100 * water[i][j][k], 1);
			glBegin(GL_QUADS);
			glVertex2d(X1, Y1);//�����`�̍����̊p
			glVertex2d(X2, Y2);//�E��
			glVertex2d(X3, Y3);//�E��
			glVertex2d(X4, Y4);//����
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//�c���C��
			X1 = converting_rate * warp[i][j][X]; Y1 = converting_rate * warp[i][j][Y];
			X2 = converting_rate * (warp[i][j][X] + yarn[i][j][X]); Y2 = converting_rate * warp[i][j][Y];
			X3 = converting_rate * (warp[i][j][X] + yarn[i][j][X]); Y3 = converting_rate * (warp[i][j][Y] + 0.5 * gap[i][j - 1][Y] + yarn[i][j][Y] + 0.5 * gap[i][j + 1][Y]);
			X4 = converting_rate * warp[i][j][X]; Y4 = converting_rate * (warp[i][j][Y] + 0.5 * gap[i][j - 1][Y] + yarn[i][j][Y] + 0.5 * gap[i][j + 1][Y]);
			/*X1 = floor(X1); Y1 = floor(Y1);
			X2 = floor(X2); Y2 = floor(Y2);
			X3 = floor(X3); Y3 = floor(Y3);
			X4 = floor(X4); Y4 = floor(Y4);*/
			glPushMatrix();
			glColor3d(1 - dyeDraw[i][j][k], 1 - dyeDraw[i][j][k], 1);
			//glColor3d(1 - c[i][j][k], 1 - c[i][j][k], 1);
			//glColor3d(1 - 100 * water[i][j][k], 1 - 100 * water[i][j][k], 1);
			glBegin(GL_QUADS);
			glVertex2d(X1, Y1);//�����`�̍����̊p
			glVertex2d(X2, Y2);//�E��
			glVertex2d(X3, Y3);//�E��
			glVertex2d(X4, Y4);//����
			glEnd();
			glPopMatrix();
		}
	}
	else if (k == 1) {//��w��
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//�c���C��
			X1 = converting_rate * warp[i][j][X]; Y1 = converting_rate * warp[i][j][Y];
			X2 = converting_rate * (warp[i][j][X] + yarn[i][j][X]); Y2 = converting_rate * warp[i][j][Y];
			X3 = converting_rate * (warp[i][j][X] + yarn[i][j][X]); Y3 = converting_rate * (warp[i][j][Y] + 0.5 * gap[i][j - 1][Y] + yarn[i][j][Y] + 0.5 * gap[i][j + 1][Y]);
			X4 = converting_rate * warp[i][j][X]; Y4 = converting_rate * (warp[i][j][Y] + 0.5 * gap[i][j - 1][Y] + yarn[i][j][Y] + 0.5 * gap[i][j + 1][Y]);
			/*X1 = floor(X1); Y1 = floor(Y1);
			X2 = floor(X2); Y2 = floor(Y2);
			X3 = floor(X3); Y3 = floor(Y3);
			X4 = floor(X4); Y4 = floor(Y4);*/
			glPushMatrix();
			glColor3d(1 - dyeDraw[i][j][k], 1 - dyeDraw[i][j][k], 1);
			//glColor3d(1 - c[i][j][k], 1 - c[i][j][k], 1);
			//glColor3d(1 - 100 * water[i][j][k], 1 - 100 * water[i][j][k], 1);
			glBegin(GL_QUADS);
			glVertex2d(X1, Y1);//�����`�̍����̊p
			glVertex2d(X2, Y2);//�E��
			glVertex2d(X3, Y3);//�E��
			glVertex2d(X4, Y4);//����
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//�����C��
			X1 = converting_rate * weft[i][j][X]; Y1 = converting_rate * weft[i][j][Y];
			X2 = converting_rate * (weft[i][j][X] + 0.5 * gap[i - 1][j][X] + yarn[i][j][X] + 0.5 * gap[i + 1][j][X]); Y2 = converting_rate * weft[i][j][Y];
			X3 = converting_rate * (weft[i][j][X] + 0.5 * gap[i - 1][j][X] + yarn[i][j][X] + 0.5 * gap[i + 1][j][X]); Y3 = converting_rate * (weft[i][j][Y] + yarn[i][j][Y]);
			X4 = converting_rate * weft[i][j][X]; Y4 = converting_rate * (weft[i][j][Y] + yarn[i][j][Y]);
			/*X1 = floor(X1); Y1 = floor(Y1);
			X2 = floor(X2); Y2 = floor(Y2);
			X3 = floor(X3); Y3 = floor(Y3);
			X4 = floor(X4); Y4 = floor(Y4);*/
			glPushMatrix();
			glColor3d(1 - dyeDraw[i][j][k], 1 - dyeDraw[i][j][k], 1);
			//glColor3d(1 - c[i][j][k], 1 - c[i][j][k], 1);
			//glColor3d(1 - 100 * water[i][j][k], 1 - 100 * water[i][j][k], 1);
			glBegin(GL_QUADS);
			glVertex2d(X1, Y1);//�����`�̍����̊p
			glVertex2d(X2, Y2);//�E��
			glVertex2d(X3, Y3);//�E��
			glVertex2d(X4, Y4);//����
			glEnd();
			glPopMatrix();
		}
	}

}