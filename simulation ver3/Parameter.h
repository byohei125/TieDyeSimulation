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
//�@�@Parameter
//�@�@Desc : ���̓p�����[�^
//---------------------------------------------------------------------------------------------------
void Parameter() {

	int cloth_setting, position_yuragi_setting, thickness_yuragi_setting;
	double densitymin, densitymax;//0<2r<m�C���Ԃ̑傫�����������̑�����菬����

	cout << "�z�̏����l��ݒ肵�܂����H�i0�F����C1�F�ݒ�1�C2�F�ݒ�2...�j"; cin >> cloth_setting;
	switch (cloth_setting) {
	case 0:
		cout << "<�J�X�^���ݒ�>" << endl;
		cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = "; cin >> yarny;
		cout << "�@����(yarnx) = "; cin >> yarnx;
		densitymin = (1 + (25.4 / yarny)) / 2; densitymax = 25.4 / yarny;
		cout << "���x [�{/inch] �F" << endl << "�@�c�����x(" << densitymin << "�`" << densitymax << ") = "; cin >> densityy;
		densitymin = ((25.4 / yarnx) + 1) / 2; densitymax = 25.4 / yarnx;
		cout << "�@�������x(" << densitymin << "�`" << densitymax << ") = "; cin >> densityx;
		cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = "; cin >> height_real;
		cout << "�@�� = "; cin >> width_real;
		cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : "; cin >> converting_rate;
		break;
	case 1:
		cout << "<�f�t�H���g�ݒ�1>" << endl;
		yarnx = 0.15; yarny = 0.15; densityx = 120; densityy = 120;
		width_real = height_real = 20;
		converting_rate = 50;
		cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
		cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
		cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
		cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : " << converting_rate << endl;
		break;
	case 2:
		cout << "<�f�t�H���g�ݒ�2>" << endl;
		yarnx = yarny = 0.15; densityx = densityy = 120;
		width_real = height_real = 100;
		converting_rate = 10;
		cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
		cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
		cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
		cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : " << converting_rate << endl;
		break;
	case 3:
		cout << "<�f�t�H���g�ݒ�3>" << endl;
		yarnx = 0.15; yarny = 0.15; densityx = 120; densityy = 120;
		width_real = height_real = 5;
		converting_rate = 200;
		cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
		cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
		cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
		cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : " << converting_rate << endl;
		break;
	}

	cout << endl << "�ʒu�̂�炬��ݒ肵�܂����H�i0�F�Ȃ��C1�F�ݒ�1...�C9�F���́j"; cin >> position_yuragi_setting;
	switch (position_yuragi_setting) {
	case 0:
		cout << "��炬�Ȃ�" << endl;
		tajudox = nx = ex = yuragix_range = 0;
		break;
	case 1:
		cout << "<�f�t�H���g�ݒ�1>" << endl;
		ax = M_PI, bx = M_PI / 7, cx = M_PI / 7;
		tajudox = 2; nx = M_PI / 2; ex = 1;
		yuragix_range = 0.1;
		cout << "�@ax = " << ax << endl;
		cout << "�@bx = " << bx << endl;
		cout << "�@cx = " << cx << endl;
		cout << "�@���d�x(tajudox) = " << tajudox << endl;
		cout << "�@���g���W��(nx) = " << nx << endl;
		cout << "�@��(ex) = " << ex << endl;
		cout << "�@��炬��(yuragix_range) = " << yuragix_range << "(" << 1 - yuragix_range << "�`" << 1 + yuragix_range << ")" << endl;
		break;	
	case 9:
		cout << "<�J�X�^���ݒ�>" << endl;
		cout << "�@���d�x(tajudox) = "; cin >> tajudox;
		cout << "�@���g���W��(nx) = "; cin >> nx;
		cout << "�@��(ex) = "; cin >> ex;
		cout << "�@��炬��(yuragix_range) = "; cin >> yuragix_range;
		cout << "�@" << 1 - yuragix_range << "�`" << 1 + yuragix_range << endl;
		break;
	}

	cout << endl << "�����̂�炬��ݒ肵�܂����H�i0�F�Ȃ��C1�F�ݒ�1...�C9�F���́j"; cin >> thickness_yuragi_setting;
	switch (thickness_yuragi_setting) {
	case 0:
		cout << "��炬�Ȃ�" << endl;
		tajudow = nw = ew = yuragiw_range = 0;
		break;
	case 1:
		cout << "<�f�t�H���g�ݒ�1>" << endl;
		tajudow = 4; nw = M_PI / 3; ew = 1;
		aw = 1.0, bw = M_PI / 3, cw = M_PI / 3;
		yuragiw_range = 0.1;
		cout << "�@aw = " << aw << endl;
		cout << "�@bw = " << bw << endl;
		cout << "�@cw = " << cw << endl;
		cout << "�@���d�x(tajudow) = " << tajudow << endl;
		cout << "�@���g���W��(nw) = " << nw << endl;
		cout << "�@��(ew) = " << ew << endl;
		cout << "�@��炬��(yuragiw_range) = " << yuragiw_range << "(" << 1 - yuragiw_range << "�`" << 1 + yuragiw_range << ")" << endl;
		break;
	case 9:
		cout << "<�J�X�^���ݒ�>" << endl;
		cout << "�@���d�x(tajudow) = "; cin >> tajudow;
		cout << "�@���g���W��(nw) = "; cin >> nw;
		cout << "�@��(ew) = "; cin >> ew;
		cout << "�@��炬��(yuragiw_range) = "; cin >> yuragiw_range;
		cout << "�@" << 1 - yuragiw_range << "�`" << 1 + yuragiw_range << endl;
		break;
	}

	cout << endl << "�z��܂�H�i0�F�܂�Ȃ��C1�F�܂�j"; cin >> fold_setting;
	if (fold_setting == 0) cout << "�܂�Ȃ�" << endl;
	else if (fold_setting == 1) {
		fold_times = 1;
		cout << "�܂�񐔂́H"; /*cin >> fold_times;*/ cout << fold_times << endl;
		for (int n = 0; n < fold_times; n++) {
			cout << "�܂���̐ݒ�F" << n << "�{��" << endl;
			//�ȉ����[��
			if (n == 0) {
				line_a[n] = 1;//-1.5
				line_b[n] = 0;//1
				line_c[n] = -94;//50
				line_a2[n] = pow(line_a[n], 2);
				line_b2[n] = pow(line_b[n], 2);
			}
			else if (n == 1) {
				line_a[n] = 0;
				line_b[n] = 1;
				line_c[n] = -94;
				line_a2[n] = pow(line_a[n], 2);
				line_b2[n] = pow(line_b[n], 2);
			}
			if (line_a[n] != 0 && line_b[n] != 0 && line_c[n] != 0) cout << line_a[n] << " * x" << line_b[n] << " * y" << " + " << line_c[n] << endl;
			else if (line_a[n] == 0) cout << line_b[n] << " * y = " << -line_c[n] << endl;
			else if (line_b[n] == 0) cout << line_a[n] << " * x = " << -line_c[n] << endl;
			else if (line_c[n] == 0) cout << line_a[n] << " * x" << line_b[n] << " * y = 0" << endl;
		}
	}

	cout << endl << "�z�͊���or�����H�i0�F�����C1�F�����j"; cin >> wet_or_dry;
	if (wet_or_dry == 0) cout << "�z�͊���" << endl;
	else cout << "�z�͎���" << endl;

	cout << endl << "�g�U�W����ݒ肵�Ă��������@"; /*cin >> D;*/
	Dy0 = 0.1;
	Dy1 = 0.1;
	Dg0 = 0.1;
	Dg1 = 0.1;
	Dg2 = 0.0001;
	D_pequal = 0.5;
	D_match = 0.3;
	cout << endl;
	cout << "Dy0 = " << Dy0 << "�@�@���� �� ������" << endl;
	cout << "Dy1 = " << Dy1 << "�@�@���� �� ������" << endl;
	cout << "Dg0 = " << Dg0 << "�@�@���� �� �������� (�g���Ȃ�)" << endl;
	cout << "Dg1 = " << Dg1 << "�@�@���� �� ��������" << endl;
	cout << "Dg2 = " << Dg2 << "�@�@�� �� ��������" << endl;
	cout << "D_pequal = " << D_pequal << "�@�@�t��" << endl;
	cout << "D_match = " << D_match << "�@�@�܂荇�킳�ꂽ�z��" << endl;

	cout << endl << "�����n�t�̉��x��ݒ肵�Ă��������i0, 20, 40, 60, 80, 100���j"; /*cin >> Temperature;*/
	Temperature = 60; cout << Temperature << "��" << endl;
	SurfaceTension = 0.06618;//60���̂Ƃ�
	//Viscosity = 0.0004658;//60���̂Ƃ�
	Viscosity = 0.04658;//60���̂Ƃ��C����͂����̃f�[�^
	/*switch (Temperature) {
	case 0:
		SurfaceTension = 0.07564; Viscosity = 0.0017917;
		break;
	case 20:
		SurfaceTension = 0.07275; Viscosity = 0.0010023;
		break;
	case 40:
		SurfaceTension = 0.06959; Viscosity = 0.0006522;
		break;
	case 60:
		SurfaceTension = 0.06618; Viscosity = 0.0004658;
		break;
	case 80:
		SurfaceTension = 0.06261; Viscosity = 0.0003550;
		break;
	case 100:
		SurfaceTension = 0.05885; Viscosity = 0.0002824;
		break;
	}*/
	
	loop_times = 100000;

}