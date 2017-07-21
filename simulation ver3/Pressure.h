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
//　　Pressure
//　　Desc : 圧迫
//---------------------------------------------------------------------------------------------------
void Pressure(int i, int j) {
	
	double var = 0.4;//ガウス関数の分散，0.64
	//まず，圧迫されているところから，圧力をひろげてく
	double rp = 10, dp;//一つの圧迫位置から影響力のある範囲の半径と，圧迫位置からの距離
	double Ap = 5;//ガウス関数の分母の係数，普通のガウス関数は2

	for (int x = i - rp; x <= i + rp; x++) {
		for (int y = j - rp; y <= j + rp; y++) {
			dp = 0;
			if (x >= 0 && y >= 0) dp = sqrt((x - i) * (x - i) + (y - j) * (y - j));//x,y<0のときは配列の用意無し
			if (dp < rp && dp > 0) {	
				p[x][y][X] += (1 / (2 * M_PI * var)) * exp(-(dp * dp) / (Ap * var));
			}
		}
	}

}

//---------------------------------------------------------------------------------------------------
//　　PressureCheck
//　　Desc : 
//---------------------------------------------------------------------------------------------------
void PressureCheck(int i, int j, int k) {

	I = i / 2; J = j / 2;	
	//布の構造の隙間（縦糸と横糸の格子の隙間）
	double gapcolor = 0;
	if (i % 2 == 1 && j % 2 == 1) {
		//周囲の糸の平均の色にする，4
		gapcolor = (p[i + 1][j + 1][k] + p[i - 1][j - 1][k] + p[i + 1][j - 1][k] + p[i - 1][j + 1][k]) / 4;
		glPushMatrix();
		glColor3d(2 - gapcolor, 2 - gapcolor, 2 - gapcolor);
		//glColor3d(1.5 - gapcolor, 1.5 - gapcolor, 1.5 - gapcolor);
		glBegin(GL_QUADS);//(I,J)座標において，1つの糸セルにつき，その右上に隙間を描画
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2));//左下
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy + GAPy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2));//右下
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy + GAPy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2) + GAPx);//右上
		glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx) + YARNx + (GAPx / 2) + GAPx);//左上
		glEnd();
		glPopMatrix();
	}

	//描画の順番で縦糸と横糸の重なりを表現するしかない
	//glVertex3dでzの値を変化させて差をつけてもうまくいかない
	if (k == 0) {//一層目
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//横糸，下
			glPushMatrix();
			glColor3d(2 - p[i][j][k], 2 - p[i][j][k], 2 - p[i][j][k]);
			//glColor3d(1.5 - p[i][j][k], 1.5 - p[i][j][k], 1.5 - p[i][j][k]);
			glBegin(GL_QUADS);
			glVertex2d((I - 1) * (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx));//長方形の左下の角
			glVertex2d((I - 1) * (YARNy + GAPy) + (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx));//右下
			glVertex2d((I - 1) * (YARNy + GAPy) + (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx) + YARNx);//右上
			glVertex2d((I - 1) * (YARNy + GAPy), (GAPx / 2) + (J - 1) * (YARNx + GAPx) + YARNx);//左上
			glEnd();
			glPopMatrix();
		}
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//縦糸，下
			glPushMatrix();
			glColor3d(2 - p[i][j][k], 2 - p[i][j][k], 2 - p[i][j][k]);
			//glColor3d(1.5 - p[i][j][k], 1.5 - p[i][j][k], 1.5 - p[i][j][k]);
			glBegin(GL_QUADS);
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy), (J - 1) * (YARNx + GAPx));//左下
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx));//右下
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy) + YARNy, (J - 1) * (YARNx + GAPx) + (YARNx + GAPx));//右上
			glVertex2d((GAPy / 2) + (I - 1) * (YARNy + GAPy), (J - 1) * (YARNx + GAPx) + (YARNx + GAPx));//左上
			glEnd();
			glPopMatrix();
		}
	}
	else if (k == 1) {//二層目
		if ((i % 4 == 0 && j % 4 == 0) || (i % 4 == 2 && j % 4 == 2)) {//縦糸，上
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
		else if ((i % 4 == 0 && j % 4 == 2) || (i % 4 == 2 && j % 4 == 0)) {//横糸，上
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