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
void Parameter(int parameter_setting) {

	int cloth_setting, yuragi_setting;
	double densitymin, densitymax;//0<2r<m�C���Ԃ̑傫�����������̑�����菬����

	switch (parameter_setting) {
	case 0://�A���P�[�g�p
		int D_enquete, Viscosity_enquete, yuragi_range_enquete;
		cout << "�A���P�[�g�p�ݒ�" << endl;
		cout << "�@�g�U�W��(1:0.0001, 2:0.001, 3:0.01)�F"; cin >> D_enquete;
		switch (D_enquete) {
		case 1:
			D = 0.0001;
			break;
		case 2:
			D = 0.001;
			break;
		case 3:
			D = 0.01;
			break;
		}
		cout << "�@�S���W��(1:0.1, 2:0.2, 3:0.5)�F"; cin >> Viscosity_enquete;
		switch (Viscosity_enquete) {
		case 1:
			Viscosity = 0.1;
			break;
		case 2:
			Viscosity = 0.2;
			break;
		case 3:
			Viscosity = 0.5;
			break;
		}
		cout << "�@��炬��(1:0.2750, 2:0.2775, 3:0.2800)�F"; cin >> yuragi_range_enquete;
		switch (yuragi_range_enquete) {
		case 1:
			yuragi_range = 0.2750;
			break;
		case 2:
			yuragi_range = 0.2775;
			break;
		case 3:
			yuragi_range = 0.2800;
			break;
		}

		cout << "�z�̐ݒ�" << endl;
		yarnx = yarny = 0.1; densityx = densityy = 200;
		width_real = height_real = 20;
		converting_rate = 40;
		cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
		cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
		cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
		cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : "; cin >> converting_rate;

		cout << "��炬�ݒ�" << endl;
		tajudo = 2; n = 3; e = 1;
		cout << "���d�x(tajudo) = " << tajudo << endl;
		cout << "���g���W��(n) = " << n << endl;
		cout << "��(e) = " << e << endl;
		cout << "��炬��(yuragi_range) = " << yuragi_range << "(" << 1 - yuragi_range << "�`" << 1 + yuragi_range << ")" << endl;

		cout << endl << "���x��ݒ肵�Ă��������i0, 20, 40, 60, 80, 100���j"; /*cin >> Temperature;*/
		Temperature = 60; cout << Temperature << "��" << endl;
		SurfaceTension = 0.06618;//60���̂Ƃ�
		loop_times = 600;

		break;

	case 1:
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
			cout << "<�f�t�H���g�ݒ�1(�A���P�[�g�Ɠ���)>" << endl;
			yarnx = yarny = 0.1; densityx = densityy = 200;
			width_real = height_real = 20;
			converting_rate = 40;
			cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
			cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
			cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
			cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : " << converting_rate << endl;
			break;
		case 2:
			cout << "<�f�t�H���g�ݒ�2>" << endl;
			yarnx = yarny = 0.2; densityx = densityy = 100;
			width_real = height_real = 20;
			converting_rate = 40;
			cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
			cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
			cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
			cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : " << converting_rate << endl;
			break;
		case 3:
			cout << "<�f�t�H���g�ݒ�3>" << endl;
			yarnx = 0.15; yarny = 0.15; densityx = 120; densityy = 120;
			width_real = height_real = 20;
			converting_rate = 50;
			cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
			cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
			cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
			cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : " << converting_rate << endl;
			break;
		case 4:
			cout << "<�f�t�H���g�ݒ�4>" << endl;
			yarnx = 0.15; yarny = 0.15; densityx = 120; densityy = 120;
			width_real = height_real = 25;
			converting_rate = 40;
			cout << "���̑��� [mm] �F" << endl << "�@�c��(yarny) = " << yarny << endl << "�@����(yarnx) = " << yarnx << endl;
			cout << "���x [�{/inch] �F" << endl << "�@�c�����x = " << densityy << endl << "�@�������x = " << densityx << endl;
			cout << "�z�̑傫�� [mm] �F" << endl << "�@�c = " << height_real << endl << "�@�� = " << width_real << endl;
			cout << "�V�~�����[�V����(�z�T�C�Y)�Ɖ�ʕ\���̔䗦�� 1 : " << converting_rate << endl;
			break;
		}

		cout << endl << "��炬�p�����[�^��ݒ肵�܂����H�i0�F����C1�F�ݒ�1...�C9�F��炬�Ȃ��j"; cin >> yuragi_setting;
		switch (yuragi_setting) {
		case 0:
			cout << "<�J�X�^���ݒ�>" << endl;
			cout << "���d�x(tajudo) = "; cin >> tajudo;
			cout << "���g���W��(n) = "; cin >> n;
			cout << "��(e) = "; cin >> e;
			cout << "��炬��" << endl << "�@yuragi_range = "; cin >> yuragi_range;
			cout << "�@" << 1 - yuragi_range << "�`" << 1 + yuragi_range << endl;
			break;
		case 1:
			cout << "<�f�t�H���g�ݒ�1(�A���P�[�g�Ɠ���)>" << yuragi_setting << ">" << endl;
			tajudo = 1; n = 1; e = 1;
			yuragi_range = 0.273;//�ł��������͗l�ƂȂ���E�l�F0.273
			cout << "���d�x(tajudo) = " << tajudo << endl;
			cout << "���g���W��(n) = " << n << endl;
			cout << "��(e) = " << e << endl;
			cout << "��炬��(yuragi_range) = " << yuragi_range << "(" << 1 - yuragi_range << "�`" << 1 + yuragi_range << ")" << endl;
			break;
		case 2:
			cout << "<�f�t�H���g�ݒ�2(�A���P�[�g�Ɠ���)>" << endl;
			tajudo = 2; n = 3; e = 1;
			yuragi_range = 0.2750;
			//yuragi_range = 0.2775;
			cout << "���d�x(tajudo) = " << tajudo << endl;
			cout << "���g���W��(n) = " << n << endl;
			cout << "��(e) = " << e << endl;
			cout << "��炬��(yuragi_range) = " << yuragi_range << "(" << 1 - yuragi_range << "�`" << 1 + yuragi_range << ")" << endl;
			break;
		case 9:
			cout << "��炬�Ȃ�" << endl;
			tajudo = 0; n = 0; e = 0; yuragi_range = 0;
			break;
		}

		cout << endl << "�z�͊���or�����H�i0�F�����C1�F�����j"; cin >> wet_or_dry;
		if (wet_or_dry == 0) cout << "�z�͊���" << endl;
		else cout << "�z�͎���" << endl;

		cout << endl << "�g�U�W����ݒ肵�Ă��������@"; /*cin >> D;*/
		Dw = 0.01;
		Dy0 = 0.01;
		Dy1 = 0.01;
		Dg0 = 0.01;
		Dg1 = 0.01;
		Dg2 = 0.0001;
		cout << endl;
		cout << "Dw = " << Dw << "�@�@���� �� �� (�����Z��)" << endl;
		cout << "Dy0 = " << Dy0 << "�@�@���� �� ������" << endl;
		cout << "Dy1 = " << Dy1 << "�@�@���� �� ������" << endl;
		cout << "Dg0 = " << Dg0 << "�@�@���� �� �������� (�g���Ȃ�)" << endl;
		cout << "Dg1 = " << Dg1 << "�@�@���� �� ��������" << endl;
		cout << "Dg2 = " << Dg2 << "�@�@�� �� ��������" << endl;

		cout << endl << "�����n�t�̉��x��ݒ肵�Ă��������i0, 20, 40, 60, 80, 100���j"; /*cin >> Temperature;*/
		Temperature = 60; cout << Temperature << "��" << endl;
		SurfaceTension = 0.06618;//60���̂Ƃ�
		//Viscosity = 0.0004658;//60���̂Ƃ�
		Viscosity = 0.04658;//60���̂Ƃ��C����͂����̃f�[�^
		loop_times = 100000;
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

		break;
	}
	cout << endl;
}