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

//---------------------------------------------------------------------------------------------------
//�@�@PressureCheck
//�@�@Desc : 
//---------------------------------------------------------------------------------------------------
void PressureCheck(int i, int j, int k) {

	I = i / 2; J = j / 2;	
	//�z�̍\���̌��ԁi�c���Ɖ����̊i�q�̌��ԁj
	double gapcolor = 0;
	if (i % 2 == 1 && j % 2 == 1) {
		//���͂̎��̕��ς̐F�ɂ���C4
		gapcolor = (p[i + 1][j + 1][k] + p[i - 1][j - 1][k] + p[i + 1][j - 1][k] + p[i - 1][j + 1][k]) / 4;
		glPushMatrix();
		glColor3d(2 - gapcolor, 2 - gapcolor, 2 - gapcolor);
		//glColor3d(1.5 - gapcolor, 1.5 - gapcolor, 1.5 - gapcolor);
		glBegin(GL_QUADS);//(I,J)���W�ɂ����āC1�̎��Z���ɂ��C���̉E��Ɍ��Ԃ�`��
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2));//����
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy + GAPy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2));//�E��
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy + GAPy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2) + GAPx);//�E��
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2) + GAPx);//����
		glEnd();
		glPopMatrix();
	}

	//�`��̏��Ԃŏc���Ɖ����̏d�Ȃ��\�����邵���Ȃ�
	//glVertex3d��z�̒l��ω������č������Ă����܂������Ȃ�
	if (k == 0) {//��w��
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//�����C��
			glPushMatrix();
			glColor3d(2 - p[i][j][k], 2 - p[i][j][k], 2 - p[i][j][k]);
			//glColor3d(1.5 - p[i][j][k], 1.5 - p[i][j][k], 1.5 - p[i][j][k]);
			glBegin(GL_QUADS);
			glVertex2d((I - 1) * (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx));//�����`�̍����̊p
			glVertex2d((I - 1) * (YARNy + GAPy) + (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx));//�E��
			glVertex2d((I - 1) * (YARNy + GAPy) + (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx) + YARNx);//�E��
			glVertex2d((I - 1) * (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx) + YARNx);//����
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//�c���C��
			glPushMatrix();
			glColor3d(2 - p[i][j][k], 2 - p[i][j][k], 2 - p[i][j][k]);
			//glColor3d(1.5 - p[i][j][k], 1.5 - p[i][j][k], 1.5 - p[i][j][k]);
			glBegin(GL_QUADS);
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy), (J - 1) * (YARNx + GAPx));//����
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx));//�E��
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx) + (YARNx + GAPx));//�E��
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy), (J - 1) * (YARNx + GAPx) + (YARNx + GAPx));//����
			glEnd();
			glPopMatrix();
		}
	}
	else if (k == 1) {//��w��
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//�c���C��
			glPushMatrix();
			glColor3d(2 - p[i][j][k], 2 - p[i][j][k], 2 - p[i][j][k]);
			//glColor3d(1.5 - p[i][j][k], 1.5 - p[i][j][k], 1.5 - p[i][j][k]);
			glBegin(GL_QUADS);
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy), (J - 1) * (YARNx + GAPx));
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx));
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx) + (YARNx + GAPx));
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy), (J - 1) * (YARNx + GAPx) + (YARNx + GAPx));
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//�����C��
			glPushMatrix();
			glColor3d(2 - p[i][j][k], 2 - p[i][j][k], 2 - p[i][j][k]);
			//glColor3d(1.5 - p[i][j][k], 1.5 - p[i][j][k], 1.5 - p[i][j][k]);
			glBegin(GL_QUADS);
			glVertex2d((I - 1) * (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx));
			glVertex2d((I - 1) * (YARNy + GAPy) + (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx));
			glVertex2d((I - 1) * (YARNy + GAPy) + (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx) + YARNx);
			glVertex2d((I - 1) * (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx) + YARNx);
			glEnd();
			glPopMatrix();
		}
	}

}