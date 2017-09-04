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
//�@�@clear
//�@�@Desc : �ϐ��̏�����
//---------------------------------------------------------------------------------------------------
void clear(void) {

	for (int i = 0; i <= NN + 1; i++) {
		for (int j = 0; j <= NN + 1; j++) {
			for (int k = 0; k <= 1; k++) {

				c[i][j][k] = 0;
				water[i][j][k] = 0;
				dye[i][j][k] = 0;
				p[i][j][k] = P_AT;//��C��1 atm
				weft[i][j][k] = 0; warp[i][j][k] = 0;
				Gap[i][j][k] = 0;
				for (int x = 0; x < 5; x++) AdjacentCellStatus[i][j][k][x] = 0;
				dCount[i][j][k] = 0;

			}
		}
	}

	cout << endl << "�ϐ��̏�����(clear)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//�@�@Cloth
//�@�@Desc : �z�̏����l�̐ݒ�
//---------------------------------------------------------------------------------------------------
void Cloth(void) {

	///////////////////////////////////////�z�̏����l�̐ݒ�///////////////////////////////////////
	gapx = (25.4 - densityx * yarnx) / (densityx - 1);//�������Ԃ̊
	gapy = (25.4 - densityy * yarny) / (densityy - 1);//�c�����Ԃ̊
	if (gapx > yarnx || gapy > yarny) cout << "�����Ӂ����̑�����茄�Ԃ̕����傫��" << endl;
	if (gapx < 0 || gapy < 0) {
		cout << "���Ԃ̑傫��gapx�܂���gapy�����̒l�̂���...�����I��" << endl;
		exit(1);
	}
	//�z���̊e���̖{���v�Z�C�ǂ������قڂ�����
	//Nx = height_real / (yarnx + gapx); Ny = width_real / (yarny + gapy);//�Ǝ��̂��傢�Ȃ��Ȍv�Z
	Nx = height_real * (densityx / 25.4);	Ny = width_real * (densityy / 25.4);//���x��p���Čv�Z
	Nx = floor(Nx); Ny = floor(Ny);
	if ((2 * Nx + 1) > NN || (2 * Ny + 1) > NN) {
		cout << "�z���̎��̖{�����ݒ肵�Ă���z��T�C�Y�������邽��...�����I��" << endl;
		exit(1);
	}
	//�`��p�̕ϊ��i�g��j
	YARNx = converting_rate * yarnx; YARNy = converting_rate * yarny;
	GAPx = converting_rate * gapx; GAPy = converting_rate * gapy;
	width = converting_rate * width_real; height = converting_rate * height_real;
	//�����_�ȉ��؎̂�
	width = floor(width); height = floor(height);
	//�l�̌ܓ��̕��@
	//Rx = floor(Rx + 0.5); Ry = floor(Ry + 0.5);

	//���Z���ƌ��ԃZ���̈�ӂ̒������`
	for (int i = 1; i <= (2 * Ny + 1); i++) {
		for (int j = 1; j <= (2 * Nx + 1); j++) {
			if (i % 2 == 0 && j % 2 == 1) {//�����Ԃ̌���
				gap[i][j][X] = yarny;//�Z����
				gap[i][j][Y] = gapx;//�Z�����s�C���Ԃ̒��a
				gap[i][j][Z] = (yarnx + yarny) / 2;//�Z������
				yarn[i][j][X] = yarny;
				yarn[i][j][Y] = 0;
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 0) {//�c���Ԃ̌���
				gap[i][j][X] = gapy;//�Z�����C���Ԃ̔��a
				gap[i][j][Y] = yarnx;//�Z�����s
				gap[i][j][Z] = (yarnx + yarny) / 2;//�Z������
				yarn[i][j][X] = 0;
				yarn[i][j][Y] = yarnx;
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 1) {//����
				gap[i][j][X] = gapy;//�Z����
				gap[i][j][Y] = gapx;//�Z�����s
				gap[i][j][Z] = (yarnx + yarny) / 2;//�Z������
				yarn[i][j][X] = 0;
				yarn[i][j][Y] = 0;
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 0 && j % 2 == 0) {//���C�z��
				gap[i][j][X] = 0;
				gap[i][j][Y] = 0;
				yarn[i][j][X] = yarny;//�Z����
				yarn[i][j][Y] = yarnx;//�Z�����s
				yarn[i][j][Z] = (yarnx + yarny) / 2;//�Z������
				capacity[i][j] = yarn[i][j][X] * yarn[i][j][Y] * yarn[i][j][Z];
			}
		}
	}

	cout << endl << "�z�̏����l�ݒ�(Cloth)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//�@�@Shibori
//�@�@Desc : �z�̈����̐ݒ�C�䂭�䂭�͍i��̐ݒ�
//---------------------------------------------------------------------------------------------------
void Shibori(void) {

	//���������ʒu�̐ݒ�
	double A1 = 0, A2 = 0, A3 = 0, A4 = 0;
	int n = 0;
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			A1 = 0.9 * 2 * Ny - i; A2 = 1.1 * 2 * Ny - i;
			A3 = 0.9 * 2 * Ny - i; A4 = 1.1 * 2 * Ny - i;
			if (j > A1 && j < A2) {
				p[i][j][X] = p[i][j][Y] = P_MAX;
				p0_i[n] = i; p0_j[n] = j;//���������ʒu���L�^
				n++;
			}
		}
	}

	//�K�E�V�A�����z
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			//��ŋL�^���������ʒu�̂Ƃ���̂݁C���͕��U���K�E�X�֐��Ōv�Z
			for (int l = 0; l <= n; l++) {
				if (p0_i[l] == i && p0_j[l] == j) {
					Pressure(i, j);
				}
				p[i][j][Y] = p[i][j][X];//�c���w�Ɖ����w�œ������͕��z
			}
		}
	}

	//���͂����ƂɎ��E���Ԃ̍�����ω�������
	double ap;//���͂ɂ��Z���̖c�����C�Z���̕��E���s�����ǂꂭ�炢�傫���Ȃ邩
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

	//��炬�ƈ��͂Ɋ�Â��ăZ���̋��e�ʂ��Čv�Z
	for (int i = 0; i <= (2 * Ny + 1); i++) {
		for (int j = 0; j <= (2 * Nx + 1); j++) {
			if (i % 2 == 0 && j % 2 == 1) {//�����Ԃ̌���
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 0) {//�c���Ԃ̌���
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 1 && j % 2 == 1) {//����
				capacity[i][j] = gap[i][j][X] * gap[i][j][Y] * gap[i][j][Z];
			}
			else if (i % 2 == 0 && j % 2 == 0) {//���C�z��
				capacity[i][j] = yarn[i][j][X] * yarn[i][j][Y] * yarn[i][j][Z];
			}
		}
	}

	cout << endl << "�i��ݒ�(Shibori)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//�@�@Dye
//�@�@Desc : �����̏����l�̐ݒ�C�����̗^����ʒu�Ȃ�
//---------------------------------------------------------------------------------------------------
void Dye(void) {
	//////////////////////////////////////�����̏����l�̐ݒ�///////////////////////////////////////
	double rA;
	//���F�O�̕z�̏�Ԃ�ύX
	if (wet_or_dry == 1) {
		for (int i = 1; i <= 2 * Ny + 1; i++) {
			for (int j = 1; j <= 2 * Nx + 1; j++) {
				for (int k = 0; k <= 1; k++) {
					water[i][j][k] = capacity[i][j] * 1.0;
				}
			}
		}
	}
	//�Ŏ΂߂ɖh�������ʒu�����E��ɉ~��ɓH��
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
	//���S�ɓH��
	/*rA = 60;
	cout << "�H���͈͂̔��a rA = " << rA << endl;
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
	cout << "�����n�t�ʁF" << dyewater << endl;
	//���S1�_�݂̂ɓH��
	/*for (int k = 0; k <= 1; k++) {
		dye[Ny + 1][Nx + 1][k] += 0.1;
		water[Ny + 1][Nx + 1][k] += 0.1;
		c[Ny + 1][Nx + 1][k] = dye[Ny + 1][Nx + 1][k] / (water[Ny + 1][Nx + 1][k] + dye[Ny + 1][Nx + 1][k]);
	}*/

	cout << "�����̏����l�ݒ�(Dye)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//�@�@ClothDraw
//�@�@Desc : �`��̂Ƃ��̃Z���̑傫���Ƃ��z�̑傫���̐ݒ�
//---------------------------------------------------------------------------------------------------
void ClothDraw(void) {

	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			if (i % 2 == 0 && j % 2 == 0) {
				weft[i][j][X] = weft[i - 2][j][X] + 0.5 * gap[i - 3][j][X] + yarn[i - 2][j][X] + 0.5 * gap[i - 1][j][X];
				weft[i][j][Y] = weft[i][j - 2][Y] + yarn[i][j - 2][Y] + gap[i][j - 1][Y];
				warp[i][j][X] = warp[i - 2][j][X] + yarn[i - 2][j][X] + gap[i - 1][j][X];
				warp[i][j][Y] = warp[i][j - 2][Y] + 0.5 * gap[i][j - 3][Y] + yarn[i][j - 2][Y] + 0.5 * gap[i][j - 1][Y];
				//�ō��[�ƍŉ��[�̐ݒ�
				weft[2][j][X] = 0;
				weft[i][2][Y] = 0.5 * gap[1][1][Y];
				warp[2][j][X] = 0.5 * gap[1][1][X];
				warp[i][2][Y] = 0;
			}
			else if (i % 2 == 1 && j % 2 == 1) {
				Gap[i][j][X] = Gap[i - 2][j][X] + gap[i - 2][j][X] + gap[i - 1][j][X];
				Gap[i][j][Y] = Gap[i][j - 2][Y] + gap[i][j - 2][Y] + gap[i][j - 1][Y];
				//�ō��[�ƍŉ��[�̐ݒ�
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

	cout << "�`��̈ʒu�ݒ�(ClothDraw)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//�@�@Yuragi
//�@�@Desc : ��炬�̕t�^
//---------------------------------------------------------------------------------------------------
void Yuragi(void) {
	//��炬�̓���
	double rmin = NN;
	double mmax = 0;
	//��炬�̕��iyuragiw�̍Œ�l�`�ő�l�j��ݒ�
	yuragixmax = 0;
	yuragiwmax = 0;
	for (int i = 1; i <= 2 * Ny + 1; i++) {//���ׂẴZ���ł�炬�v�Z�Ci = 0�ɂ����Ă�炬��0
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

	for (int i = 1; i <= 2 * Ny + 1; i++) {//i = 0�ɂ����Čv�Z����Ƃ�炬��0�ɂȂ�
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
	
	cout << "��炬�ݒ�(Yuragi)...COMPLETE" << endl;

}

//---------------------------------------------------------------------------------------------------
//�@�@Display1
//�@�@Desc : �E�B���h�E1�ւ̕`��
//---------------------------------------------------------------------------------------------------
void Display1(void) {
	auto startTime = std::chrono::system_clock::now();

	//��ʃN���A
	glClear(GL_COLOR_BUFFER_BIT);

	//���F�v�Z
	double dyemax = 0;
	cout << endl << "���F�v�Z��";
	cout << endl << "Press 'q' to get out of the loop" << endl << "Press 't' to get interim output" << endl;
	for (t = 1; t <= loop_times; t++) {
		dTerminator = 0;
		bcTerminator = 0;

		//���ڃZ���̋ߖT�̃Z���̏��(���� or ����)���m�F
		for (int i = 1; i <= 2 * Ny + 1; i++) {
			for (int j = 1; j <= 2 * Nx + 1; j++) {
				for (int k = 0; k <= 1; k++) {
					//�אڃZ���̊����E�����𔻕�
					if (water[i + 1][j][k] >= capacity[i + 1][j] * 1.0) AdjacentCellStatus[i][j][k][0] = 1;
					if (water[i - 1][j][k] >= capacity[i - 1][j] * 1.0) AdjacentCellStatus[i][j][k][1] = 1;
					if (water[i][j + 1][k] >= capacity[i][j + 1] * 1.0) AdjacentCellStatus[i][j][k][2] = 1;
					if (water[i][j - 1][k] >= capacity[i][j - 1] * 1.0) AdjacentCellStatus[i][j][k][3] = 1;
					//�אڃZ������������Ă��邩
					if (p[i + 1][j][k] >= P_MAX) AdjacentCellStatus[i][j][k][0] = 2;
					if (p[i - 1][j][k] >= P_MAX) AdjacentCellStatus[i][j][k][1] = 2;
					if (p[i][j + 1][k] >= P_MAX) AdjacentCellStatus[i][j][k][2] = 2;
					if (p[i][j - 1][k] >= P_MAX) AdjacentCellStatus[i][j][k][3] = 2;
					//���̑w�̏�Ԃ͌��ݖ��l��(2017/06/15)
					/*if (k == 0 && water[i][j][k + 1] >= capacity[i][j] * 1.0) AdjacentCellStatus[i][j][k][4] = 1;
					else if (k == 1 && water[i][j][k - 1] >= capacity[i][j] * 1.0) AdjacentCellStatus[i][j][k][4] = 1;*/
					dCount[i][j][k] = 0;
					bcCount[i][j][k] = 0;
					d_pequalCount[i][j][k] = 0;
				}
			}
		}

		if (d_pequalTerminator > 0) {//t = t-1�̂Ƃ��̌��ʂ���
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
					if (dCount[i][j][k] != 0) dTerminator++;//���ׂẴZ���Ŋg�U�v�Z����Ȃ����dTerminator = 0
					if (bcCount[i][j][k] != 0) bcTerminator++;//���ׂẴZ����Buras�̎��E�э׊Ǎ�p�̌v�Z������Ȃ����bcTerminator = 0
					if (d_pequalCount[i][j][k] != 0) d_pequalTerminator++;
				}
			}
		}
		if (t % 100 == 0) {
			auto endTime = std::chrono::system_clock::now();
			auto timeSpan = endTime - startTime;
			cout << "t = " << t << " (" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms])" << 
				endl << "�@dTerminator = " << dTerminator << ", bcTerminator = " << bcTerminator << ", d_pequalTerminator = " << d_pequalTerminator << endl;
		}
		if (dTerminator == 0 && bcTerminator == 0 && d_pequalTerminator == 0) {
			cout << "t = " << t << endl;
			break;
		}

		//����̃L�[���������ꂽ�Ƃ��̏���
		char buf;
		if (_kbhit()) {
			buf = _getch();
			//q�L�[���������ꂽ�Ƃ��Ƀ��[�v�����o��
			if (buf == 'q') {
				cout << "t = " << t << "�ŋ����I��" << endl;
				break;
			}
			//t�L�[���������ꂽ�Ƃ��ɓr�����ʉ摜�o��
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
				WriteBitmap("Simulation output(�r���o��).bmp", width, height);
				auto endTime = std::chrono::system_clock::now();
				auto timeSpan = endTime - startTime;
				cout << "�r���o�߉摜�o��: t = " << t << " (" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms])" << endl;
			}
		}

	}

	cout << "���F...COMPLETE" << endl;

	//�����ʌv�Z
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
				if (capacity[i][j] != 0) {
					dyeDraw[i][j][k] = dye[i][j][k] / capacity[i][j] * 2;
				}
			}
		}
	}

	//�`��
	cout << endl << "�`��" << endl;
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

	cout << "�`��...COMPLETE" << endl;
	WriteBitmap("Simulation output.bmp", width, height);
	cout << "�ۑ�...COMPLETE" << endl;

	glFlush();

	auto endTime = std::chrono::system_clock::now();
	auto timeSpan = endTime - startTime;
	cout << "��������:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;

}

//---------------------------------------------------------------------------------------------------
//�@�@Display2
//�@�@Desc : �E�B���h�E2�ւ̕`��
//---------------------------------------------------------------------------------------------------
void Display2(void) {

	//��ʃN���A
	glClear(GL_COLOR_BUFFER_BIT);

	//���͕��z
	for (int i = 1; i <= 2 * Ny + 1; i++) {
		for (int j = 1; j <= 2 * Nx + 1; j++) {
			for (int k = 0; k <= 1; k++) {
			}
		}
	}
	glFlush();
}

//---------------------------------------------------------------------------------------------------
//�@�@resize
//�@�@Desc : 
//---------------------------------------------------------------------------------------------------
void resize(int w, int h) {
	/* �E�B���h�E�S�̂��r���[�|�[�g�ɂ��� */
	glViewport(0, 0, w, h);
	/* �ϊ��s��̏����� */
	glLoadIdentity();
	//�\��������W�͈͂�(0, 0, -1)~(400, 400, 1)�Ƃ���
	glOrtho(0, w, 0, h, -1.0, 1.0);
}

//----------------------------------------------------------------------------------------------------
//�@�@Initialize
//�@�@Desc : �����������C�w�i�F
//----------------------------------------------------------------------------------------------------
void Initialize(void) {
	//glClearColor(R, G, B, ���l�F�s�����x�i0 �œ���, 1 �ŕs�����j);
	glClearColor(1.0, 1.0, 1.0, 1.0);
}

//----------------------------------------------------------------------------------------------------
//�@�@main
//�@�@Desc : ���C���G���g���[�|�C���g
//----------------------------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
	
	int WinID[2];

	Parameter();//���̓p�����[�^
	clear();
	Cloth();
	ClothDraw();//��炬�Ȃ��ŕ`��ʒu���v�Z
	Shibori();
	Yuragi();	
	Dye();

	cout << "���ԁF" << endl << "�@�c���� gapy = " << gapy << endl;
	cout << "�@������ gapx = " << gapx << endl;
	cout << "�z���̎��̖{���F" << endl << "�@�c�� Ny = " << Ny << "�@(�Z���� 2Ny+1 = " << 2 * Ny + 1 << ")"
		<< endl << "�@���� Nx = " << Nx << "�@(�Z���� 2Nx+1 = " << 2 * Nx + 1 << ")" << endl;
	cout << "�Z���̕W���̐ρF" << endl << "�@���F" << yarnx * yarny * (yarnx + yarny) / 2 << endl;
	cout << "�@����(�c��)�F" << gapx * yarny * (yarnx + yarny) / 2 << endl;
	cout << "�z�̋��e��(�̐ρC��炬�Ȃ�)�F" << ((yarny + gapy) * Nx) * ((yarnx + gapx) * Ny) * (yarnx + yarny) << endl;
	/*cout << endl << "YARNx = " << YARNx << ", YARNy = " << YARNy << ", GAPx = " << GAPx << ", GAPy = " << GAPy << endl;*/
	cout << "width = " << width << ", height = " << height << endl;
	cout << "loop times = ";	cout << loop_times << endl;
	//cin >> loop_times;

	glutInit(&argc, argv);//glut�̏�����
	//�E�B���h�E1
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);//�l�̓s�N�Z��
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	WinID[0] = glutCreateWindow("Tie Dye simulation");
	glutDisplayFunc(Display1);
	glutReshapeFunc(resize);
	glutMouseFunc(mouse);
	Initialize();

	//�E�B���h�E2
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