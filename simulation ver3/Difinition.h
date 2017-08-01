#pragma once
#define NN 1000//�z�̈�ӂ̒����i��+���ԁj�̍ő�l�C2*Nx��2*Ny���傫���Ȃ��Ă͂Ȃ�Ȃ�
int Nx, Ny;//�z���̉����̖{���C�c���̖{��
int width, height;//�E�B���h�E�傫���E���W
int converting_rate;//�V�~�����[�V�����i�v�Z�C���ۂ̕z�T�C�Y�j�Ɖ�ʕ\���̕ϊ��䗦
int width_real, height_real;//���ۂ̕z�̃T�C�Y�C10 cm�Ƃ�

double c[NN + 2][NN + 2][2];//�Z���̐����̔Z�x�i0�F�܎�, 1�F�o���j�it)
//�Z�x�̊g�U�v�Z�̂Ƃ��g�p
double dcN;//�� �� ����
double dcS;
double dcW;
double dcE;
double dcN_gap;//���� �� ����
double dcS_gap;
double dcW_gap;
double dcE_gap;
double dcN_yarn;//���� �� ��
double dcS_yarn;
double dcW_yarn;
double dcE_yarn;

double dye[NN + 2][NN + 2][2];//������
double ddye;//�����̈ړ���
double ddye1;//�����̈ړ���

double water[NN + 2][NN + 2][2];//������
double dwater;//�ړ����鐅���C�Z���̖O�a�ʂ𒴂�����
double dwater1;//dwater��1�Z�����ɒ����CPlain.h�Q��

double dyewater;//�����n�t��
double ddyewater;//�Z���v�Z�ňړ���������n�t��
double ddyewater1;//�Z���v�Z�ňړ���������n�t��

//Buras�̎�
double dq;//���ԓ�����̋z����
double dq1;//���ԓ�����̋z����
double initialVelocity = 0.1;//�����z�����x
double elapsedTime;//�o�ߎ���

//�g�U�W��
double D;
//���� �� ���ւ̊g�U
double DyN;//���ڃZ���́�
double DyS;//��
double DyW;//��
double DyE;//��
double Dy0;//��������
double Dy1;//��������
//���� �� ���� or ���ւ̊g�U
double DgN;//��
double DgS;//��
double DgW;//��
double DgE;//��
double Dg0;//���� �� ���ԁC�������Ԃ�
double Dg1;//���� �� ���ԁC�������Ԃ�
double Dg2;//�� �� ���ԁC�������Ԃ�

//���� �� ���ɂ����鐅���Z���̌W��
double Dw;

//��炬
double yuragix[NN][NN];
double yuragiy[NN][NN];
double yuragiw[NN][NN];
double ax, bx, cx;
double ay, by, cy;
double aw, bw, cw;
int tajudox, tajudow;
double nx, ex, nw, ew;
double fx1, fx2, fy1, fy2, fw1, fw2;
double yuragixmax, yuragixmin;
double yuragiymax, yuragiymin;
double yuragiwmax, yuragiwmin;
double yuragix_range;//��炬�̕��C1-yuragi_range ~ 1+yuragi_range�C0.8~1.2�݂�����
double yuragiw_range;//��炬�̕��C1-yuragi_range ~ 1+yuragi_range�C0.8~1.2�݂�����
double yuragix_range_setting;//��炬�ʂ�yuragi_range�ɂ��邽��
double yuragiw_range_setting;//��炬�ʂ�yuragi_range�ɂ��邽��

//����
double p[NN + 2][NN + 2][2];//����
int p0_i[NN * 100], p0_j[NN * 100];//�����ʒu�C�ŏ��ɕz���������Ă���Ƃ���

//�э׊Ǎ�p
int Temperature;//���x
double SurfaceTension;//�\�ʒ���[N/m]
const double ContactAngle = 0;//�ڐG�p��
const double LiquidDensity = 1000;//�t�̖̂��x[kg/m^3]
double Viscosity;//�t�̂̔S�x�i�S���W���j[Pa*s]
double r;//�э׊Ǎ�p�̐Z���������v�Z����ۂ̌��Ԃ̔��a
double rx;//x(i)�����̐Z�������v�Z�̂Ƃ��̌��Ԕ��a
double ry;//y(j)�����̐Z�������v�Z�̂Ƃ��̌��Ԕ��a
double h;//�э׊Ǎ�p�ɂ��Z������
double h1;//�v�Z�p
double hx;//�э׊Ǎ�p�ɂ��x(i)�����̐Z������
double hx1;//�v�Z�p
double hy;//�э׊Ǎ�p�ɂ��y(j)�����̐Z������
double hy1;//�v�Z�p
int hpositivex, hnegativex, hpositivey, hnegativey;//�Z���Z�������J�E���g


//�z�̐ݒ�
double densityx;//�����̖��x�i130~253�C1�C���`=25.4cm���j
double densityy;//�c���̖��x
double gap[NN + 2][NN + 2][3];//���Ԃ̌��Ԃ̑傫��[mm]�C0�F�Z����x�C1�F�Z�����sy�C2�F�Z������z
double gapx;//�����Ԃ̌��Ԃ̒��a�̊�l[mm]
double gapy;//�c���Ԃ̌��Ԃ̒��a�̊�l[mm]
double yarn[NN + 2][NN + 2][3];//���̑��� = �Z���̈�ӂ̒���[mm]�C0�F�Z����x�C1�F�Z�����sy�C2�F�Z������z
double yarnx;//�����̑����̊�l�i������Ԏ��蓾���鑾���j[mm]
double yarny;//�c���̑����̊�l�i������Ԏ��蓾���鑾���j[mm]
double capacity[NN + 2][NN + 2];//�Z���̖O�a�ʁi���e�ʁC�̐ρj[mm^3]


//�v�Z�p�̕ϐ�
double dt = 0.005;//1�^�C���X�e�b�v�Ɍo�߂��鎞��0.005
double c0 = 0;//���������Z�x,0.4
double dye0 = 0.2;//����������,0.2
double water0 = 0.3;//����������,0.3
int t = 0;//�v�Z���񂵂��񐔁C����
int loop_times = 0;
const int X = 0;
const int Y = 1;
const int Z = 2;
int WetDryFlag[NN + 2][NN + 2][2][5];//�l��0�̂Ƃ������E1�̂Ƃ������C0:���ڃZ���̉E�̃Z���C1:���C2:��C3:���C4:���̑w(���Ԃɂ�����)
int dTerminator = 0;//���F�I������(�g�U�v�Z)
int bcTerminator = 0;//���F�I������(Buras�̎��Ɩэ׊Ǎ�p)
int dCount[NN + 2][NN + 2][2];//�g�U�v�Z���N���Ă��邩�m�F�p
int bcCount[NN + 2][NN + 2][2];//Buras�̎��Ɩэ׊Ǎ�p�̌v�Z���N���Ă��邩�m�F�p
int wet_or_dry;//������Ԃɂ����ĕz�͎�����������(0�F�����C1�F����)


//�`��p
int YARNx, YARNy;//�y�`��p�z���Z���̈�Ӂi���a�j
double GAPx, GAPy;//�y�`��p�z���ԃZ���̈�Ӂi���a�j
int I, J;
double weft[NN + 2][NN + 2][2];//�`��J�n���W�C�����Z���̍����p�̍��W�l�C0�FX���W�C1�FY���W
double warp[NN + 2][NN + 2][2];//�`��J�n���W�C�c���Z���̍����p�̍��W�l�C0�FX���W�C1�FY���W
double Gap[NN + 2][NN + 2][2];//�`��J�n���W�C�z�̌��ԕ����̎l�p�`�̍����p�̍��W�l�C0�FX���W�C1�FY���W
double weft1[NN + 2][NN + 2][2];//�`��J�n���W�C�����Z���̍����p�̍��W�l�C0�FX���W�C1�FY���W
double warp1[NN + 2][NN + 2][2];//�`��J�n���W�C�c���Z���̍����p�̍��W�l�C0�FX���W�C1�FY���W
double Gap1[NN + 2][NN + 2][2];//�`��J�n���W�C�z�̌��ԕ����̎l�p�`�̍����p�̍��W�l�C0�FX���W�C1�FY���W
double dyeDraw[NN + 2][NN + 2][2];//�`��p�����l

//�D���
int weavepattern = 0;//Plain = 1, Twill = 2�CSateen = 3

//---------------------------------------------------------------------------------------------------
//�@�@InitRand
//�@�@Desc : �����̏�����
//---------------------------------------------------------------------------------------------------
inline void InitRand(void)
{
	srand((unsigned int)time(NULL));
}

//---------------------------------------------------------------------------------------------------
//�@�@Rand_mod
//�@�@Desc : �����C�g���₷���悤�ɒ����Ca:�����͈́i���鐔�C���܂�j�Cb:�����̌�����
//---------------------------------------------------------------------------------------------------
inline int Rand_mod(int a, double b)
{
	return rand() % a * b;//rand()��0~RAND_MAX�܂�
}

//---------------------------------------------------------------------------------------------------
//�@�@Rand_range
//�@�@Desc : �����C�����͈́Fmin~max
//---------------------------------------------------------------------------------------------------
inline int Rand_range(double min, double max)
{
	return min + (double)(rand()*(max - min + 1.0) / (1.0 + RAND_MAX));
}

