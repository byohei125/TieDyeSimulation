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
//�@�@Pressure
//�@�@Desc : ����
//---------------------------------------------------------------------------------------------------
void Pressure(int i, int j) {
	
	double var = 0.4;//�K�E�X�֐��̕��U�C0.64
	//�܂��C��������Ă���Ƃ��납��C���͂��Ђ낰�Ă�
	double rp = 10, dp;//��̈����ʒu����e���͂̂���͈͂̔��a�ƁC�����ʒu����̋���
	double Ap = 5;//�K�E�X�֐��̕���̌W���C���ʂ̃K�E�X�֐���2

	for (int x = i - rp; x <= i + rp; x++) {
		for (int y = j - rp; y <= j + rp; y++) {
			dp = 0;
			if (x >= 0 && y >= 0) dp = sqrt((x - i) * (x - i) + (y - j) * (y - j));//x,y<0�̂Ƃ��͔z��̗p�Ӗ���
			if (dp < rp && dp > 0) {	
				p[x][y][X] += (1 / (2 * M_PI * var)) * exp(-(dp * dp) / (Ap * var));
			}
		}
	}

}
