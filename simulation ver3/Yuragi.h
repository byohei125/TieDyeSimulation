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
//　　YuragiCal
//　　Desc : ゆらぎ関数
//---------------------------------------------------------------------------------------------------
void YuragiCal(int i, int j) {
	double I, J;
	I = i / M_PI / 4;
	J = j / M_PI / 4;
	ax = M_PI, bx = M_PI / 7, cx = M_PI / 7;
	aw = 1.0, bw = M_PI / 3, cw = M_PI / 3;

	for (int l = 1; l <= tajudox; l++) {
		fx1 = sin(I + ax*sin(bx*J) + (l - 1)*cx);
		fx2 = sin(J + ax*sin(bx*I) + (l - 1)*cx);
		yuragix[i][j] += (sin(pow(nx, (l - 1)) * fx1) + sin(pow(nx, (l - 1)) * fx2)) / pow(nx, ex * (l - 1));
		yuragixmax = (yuragix[i][j] >= yuragixmax) ? yuragix[i][j] : yuragixmax;
	}
	for (int l = 1; l <= tajudow; l++) {
		fw1 = sin(I + aw*sin(bw*J) + (l - 1)*cw);
		fw2 = sin(J + aw*sin(bw*I) + (l - 1)*cw);
		yuragiw[i][j] += (sin(pow(nw, (l - 1)) * fw1) + sin(pow(nw, (l - 1)) * fw2)) / pow(nw, ew * (l - 1));
		yuragiwmax = (yuragiw[i][j] >= yuragiwmax) ? yuragiw[i][j] : yuragiwmax;
	}

}