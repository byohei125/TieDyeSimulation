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
