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
//　　Parameter
//　　Desc : 入力パラメータ
//---------------------------------------------------------------------------------------------------
void Parameter() {

	int cloth_setting, yuragi_setting;
	double densitymin, densitymax;//0<2r<m，隙間の大きさが正かつ糸の太さより小さい

	cout << "布の初期値を設定しますか？（0：する，1：設定1，2：設定2...）"; cin >> cloth_setting;
	switch (cloth_setting) {
	case 0:
		cout << "<カスタム設定>" << endl;
		cout << "糸の太さ [mm] ：" << endl << "　縦糸(yarny) = "; cin >> yarny;
		cout << "　横糸(yarnx) = "; cin >> yarnx;
		densitymin = (1 + (25.4 / yarny)) / 2; densitymax = 25.4 / yarny;
		cout << "密度 [本/inch] ：" << endl << "　縦糸密度(" << densitymin << "〜" << densitymax << ") = "; cin >> densityy;
		densitymin = ((25.4 / yarnx) + 1) / 2; densitymax = 25.4 / yarnx;
		cout << "　横糸密度(" << densitymin << "〜" << densitymax << ") = "; cin >> densityx;
		cout << "布の大きさ [mm] ：" << endl << "　縦 = "; cin >> height_real;
		cout << "　横 = "; cin >> width_real;
		cout << "シミュレーション(布サイズ)と画面表示の比率は 1 : "; cin >> converting_rate;
		break;
	case 1:
		cout << "<デフォルト設定1>" << endl;
		yarnx = 0.15; yarny = 0.15; densityx = 120; densityy = 120;
		width_real = height_real = 20;
		converting_rate = 50;
		cout << "糸の太さ [mm] ：" << endl << "　縦糸(yarny) = " << yarny << endl << "　横糸(yarnx) = " << yarnx << endl;
		cout << "密度 [本/inch] ：" << endl << "　縦糸密度 = " << densityy << endl << "　横糸密度 = " << densityx << endl;
		cout << "布の大きさ [mm] ：" << endl << "　縦 = " << height_real << endl << "　横 = " << width_real << endl;
		cout << "シミュレーション(布サイズ)と画面表示の比率は 1 : " << converting_rate << endl;
		break;
	case 2:
		cout << "<デフォルト設定2>" << endl;
		yarnx = yarny = 0.15; densityx = densityy = 120;
		width_real = height_real = 100;
		converting_rate = 10;
		cout << "糸の太さ [mm] ：" << endl << "　縦糸(yarny) = " << yarny << endl << "　横糸(yarnx) = " << yarnx << endl;
		cout << "密度 [本/inch] ：" << endl << "　縦糸密度 = " << densityy << endl << "　横糸密度 = " << densityx << endl;
		cout << "布の大きさ [mm] ：" << endl << "　縦 = " << height_real << endl << "　横 = " << width_real << endl;
		cout << "シミュレーション(布サイズ)と画面表示の比率は 1 : " << converting_rate << endl;
		break;
	case 3:
		cout << "<デフォルト設定3>" << endl;
		yarnx = 0.15; yarny = 0.15; densityx = 120; densityy = 120;
		width_real = height_real = 25;
		converting_rate = 40;
		cout << "糸の太さ [mm] ：" << endl << "　縦糸(yarny) = " << yarny << endl << "　横糸(yarnx) = " << yarnx << endl;
		cout << "密度 [本/inch] ：" << endl << "　縦糸密度 = " << densityy << endl << "　横糸密度 = " << densityx << endl;
		cout << "布の大きさ [mm] ：" << endl << "　縦 = " << height_real << endl << "　横 = " << width_real << endl;
		cout << "シミュレーション(布サイズ)と画面表示の比率は 1 : " << converting_rate << endl;
		break;
	}

	cout << endl << "ゆらぎパラメータを設定しますか？（0：する，1：設定1...，9：ゆらぎなし）"; cin >> yuragi_setting;
	switch (yuragi_setting) {
	case 0:
		cout << "<カスタム設定>" << endl;
		cout << "位置ゆらぎ" << endl;
		cout << "　多重度(tajudox) = "; cin >> tajudox;
		cout << "　周波数係数(nx) = "; cin >> nx;
		cout << "　ε(ex) = "; cin >> ex;
		cout << "　ゆらぎ幅(yuragix_range) = "; cin >> yuragix_range;
		cout << "　" << 1 - yuragix_range << "〜" << 1 + yuragix_range << endl;
		cout << "太さゆらぎ" << endl;
		cout << "　多重度(tajudow) = "; cin >> tajudow;
		cout << "　周波数係数(nw) = "; cin >> nw;
		cout << "　ε(ew) = "; cin >> ew;
		cout << "　ゆらぎ幅(yuragiw_range) = "; cin >> yuragiw_range;
		cout << "　" << 1 - yuragiw_range << "〜" << 1 + yuragiw_range << endl;
		break;
	case 1:
		cout << "<デフォルト設定1>" << endl;
		tajudox = 2; nx = M_PI / 2; ex = 1;
		yuragix_range = 0.03;
		tajudow = 4; nw = M_PI / 3; ew = 1;
		yuragiw_range = 0.03;
		cout << "位置ゆらぎ" << endl;
		cout << "　多重度(tajudox) = " << tajudox << endl;
		cout << "　周波数係数(nx) = " << nx << endl;
		cout << "　ε(ex) = " << ex << endl;
		cout << "　ゆらぎ幅(yuragix_range) = " << yuragix_range << "(" << 1 - yuragix_range << "〜" << 1 + yuragix_range << ")" << endl;
		cout << "太さゆらぎ" << endl;
		cout << "　多重度(tajudow) = " << tajudow << endl;
		cout << "　周波数係数(nw) = " << nw << endl;
		cout << "　ε(ew) = " << ew << endl;
		cout << "　ゆらぎ幅(yuragiw_range) = " << yuragiw_range << "(" << 1 - yuragiw_range << "〜" << 1 + yuragiw_range << ")" << endl;
		break;
	case 9:
		cout << "ゆらぎなし" << endl;
		tajudox = nx = ex = 0;
		tajudow = nw = ew = 0;
		break;
	}

	cout << endl << "布は乾燥or湿潤？（0：乾燥，1：湿潤）"; cin >> wet_or_dry;
	if (wet_or_dry == 0) cout << "布は乾燥" << endl;
	else cout << "布は湿潤" << endl;

	cout << endl << "拡散係数を設定してください　"; /*cin >> D;*/
	Dw = 0.1;
	Dy0 = 0.1;
	Dy1 = 0.1;
	Dg0 = 0.1;
	Dg1 = 0.1;
	Dg2 = 0.0001;
	cout << endl;
	cout << "Dw = " << Dw << "　　隙間 → 糸 (水分浸透)" << endl;
	cout << "Dy0 = " << Dy0 << "　　隙間 → 乾燥糸" << endl;
	cout << "Dy1 = " << Dy1 << "　　隙間 → 湿潤糸" << endl;
	cout << "Dg0 = " << Dg0 << "　　隙間 → 乾燥隙間 (使われない)" << endl;
	cout << "Dg1 = " << Dg1 << "　　隙間 → 湿潤隙間" << endl;
	cout << "Dg2 = " << Dg2 << "　　糸 → 湿潤隙間" << endl;

	cout << endl << "染料溶液の温度を設定してください（0, 20, 40, 60, 80, 100℃）"; /*cin >> Temperature;*/
	Temperature = 60; cout << Temperature << "℃" << endl;
	SurfaceTension = 0.06618;//60℃のとき
	//Viscosity = 0.0004658;//60℃のとき
	Viscosity = 0.04658;//60℃のとき，これはうそのデータ
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

}