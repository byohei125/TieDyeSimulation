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
//　　Yuragi
//　　Desc : ゆらぎ関数
//---------------------------------------------------------------------------------------------------
void Yuragi(int i, int j) {

	for (tajudo_i = 1; tajudo_i <= tajudo; tajudo_i++) {
		//w
		aw = 1.0, bw = 1.0, cw = 0.1;
		fw1 = sin(i + aw*sin(bw*j) + (tajudo_i - 1)*cw);
		fw2 = sin(j + aw*sin(bw*i) + (tajudo_i - 1)*cw);
		yuragiw[i][j] += (sin(pow(n, (tajudo_i - 1)) * fw1) + sin(pow(n, (tajudo_i - 1)) * fw2)) / pow(n, e * (tajudo_i - 1));

		yuragiwmax = (yuragiw[i][j] >= yuragiwmax) ? yuragiw[i][j] : yuragiwmax;//これは絶対必要!!
		//yuragiwmin = (yuragiw[i][j] <= yuragiwmin) ? yuragiw[i][j] : yuragiwmin;

	}

}