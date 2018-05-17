

#include "stdafx.h"
#include "resource.h"




#define WM_NEWDATA (WM_USER+1)
UINT thread_id;

	CToolBarCtrl Bar;
	CMenu menu;


#include "vector.h"



	PARTICLES p(1);
	ATOMS A;


	REAL scale;
	REAL Tx,Ty,Tz;
	int nsteps;
	REAL radius;
	REAL dt;
	REAL zmax;

	bool stoped, ready;


	void ComputeEuler ();
	void ComputeRK2 ();
	void ComputeRK3 ();
	void ComputeRK4 ();
	void ComputeRKCK ();
	void (*Compute) ();

	void TimerInteractive ();
	void TimerNonInteractive ();
	void (*timer) ();


	
	DWORD st,fn;
	REAL _C, C, Cek, Cpot;
	REAL E0, Ek, Ep, SummREE, e1, e2;
	const REAL iatf     = 3.88005260641761e10;
	const REAL halfatf2 = 3.3211980678119621748e-22;

	VECTOR Track[1000];
	DWORD jj, ii;



  REAL hh, h6;
  REAL t, tsav, Tsav, h;
  int nok, nbad;

  static
  REAL b[6][5] = {0,               0,           0,             0,                0,
                  0.2,             0,           0,             0,                0,
                  3.0/40.0,        9.0/40.0,    0,             0,                0,
                  0.3,             -0.9,        1.2,           0,                0,
                  -11.0/54.0,      2.5,         -70.0/27.0,    35.0/27.0,        0,
                  1631.0/55296.0 , 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0
		};
  static REAL
        c1      = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
        dc1     = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0, dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6 = c6-0.25;

  static
  REAL c[6]     = {37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0};


  static void initRK()
  {
    h = dt;
    hh = dt/2.0;
    h6 = dt/6.0;
  }

//===========================================================
//===========================================================


REAL Epot(REAL r)
{
  return (1/r - pow(3+r*r, -0.5));
}

REAL Energy()
{
  int i;
  VECTOR r;
  REAL modul_r;

  //Compute kinetic & pot. energy
  Ek = Cek*Length2(p.vel[0]);
  Ep = 0;
  for (i = 0; i <= 349; i++)
  {
    r = p.pos[0] - A.pos[i];
    modul_r = Length(r);
    if (modul_r < radius)
      Ep += Epot(modul_r);
  }
  Ep = Cpot*Ep;





  return Ep + Ek;
}

REAL EnergyExact()
{
  int i;
  VECTOR r;
  REAL modul_r;

  //Compute kinetic & pot. energy
  Ek = Cek*Length2(p.vel[0]);
  Ep = 0;
  for (i = 0; i <= 349; i++)
  {
    r = p.pos[0] - A.pos[i];
    modul_r = Length(r);
    Ep += Epot(modul_r);
  }
  Ep = Cpot*Ep;
  return Ep + Ek;
}






  REAL adt(REAL r)
{
  REAL a, b, c;

  a = r*r;
  b = 3+a;
  c = pow(b, 1.5);
  a = 1/b + 3/(a*b) - r/c;
  return -C*a;
}


REAL ah(REAL r)
{
  REAL a, b, c;

  a = r*r;
  b = 3+a;
  c = pow(b, 1.5);
  a = 1/b + 3/(a*b) - r/c;
  return -_C*a*h;
}


void Euler(int i) // where "i" is a particle index
{
  int j;
  VECTOR r;
  REAL modul_r;
  VECTOR delta_v, delta_x;


    for (j = 0; j <= 349; j++)
    {
      r = p.pos[i] - A.pos[j];
      modul_r = Length(r);
      if (modul_r < radius)
      {
        delta_v = Normalize(r)*adt(modul_r);
        p.vel[i] += delta_v;
      }
    }
    delta_x = p.vel[i]*dt;
    p.pos[i] += delta_x;
}


void ComputeEuler()
{
  int i;
  

  for (i = 1; i <= nsteps; i++)
    if (!p.finished[0])
    {
      Euler(0);

      if (p.pos[0].z > zmax)
      {
        p.finished[0] = true;
        stoped = true;
      }
    }

  if (!p.finished[0])
  {
    Track[jj] = p.pos[0];
    jj += 1;
  }
}



void rk2(int i)
{
  int j;
  VECTOR r;
  VECTOR k1, k2;
  VECTOR v, x;
  VECTOR oldpos, oldvel;
  REAL modul_r;

  oldpos = p.pos[i];
  oldvel = p.vel[i];

  Null(x);
  Null(v);

  for (j = 0; j <= 349; j++)
  {
    r = oldpos - A.pos[j];
    modul_r = Length(r);
    if (modul_r < radius)
    {
      k1 = Normalize(r)*adt(Length(r)); //First step.
      r += oldvel*hh;
      k2 = Normalize(r)*adt(Length(r));
      v += k2;
      x += k1*hh;
    }
  }
  p.vel[i] += v;
  p.pos[i] += oldvel*dt + x;
}


void ComputeRK2()
{
  int i;

  for (i = 1; i <= nsteps; i++)
    if (!p.finished[0])
    {
      rk2(0);
      if (p.pos[0].z > zmax)
      {
        p.finished[0] = true;
        stoped = true;
        /*fn = GetTickCount(); inttostr(fn-st);*/
      }
    }
  if (!p.finished[0])
  {
    Track[jj] = p.pos[0];
    jj += 1;
  }
}




void rk3(int i)
{
  int j;
  VECTOR r, r2;
  VECTOR k1, k2, k3;
  VECTOR v,x;
  VECTOR oldpos, oldvel, oldveldt;
  REAL modul_r;

  Null(x);
  Null(v);

  oldpos = p.pos[i];
  oldvel = p.vel[i];
  oldveldt = oldvel*dt;

  for (j = 0; j <= 349; j++)
  {
    r = oldpos - A.pos[j];
    modul_r = Length(r);
    if (modul_r < radius)
    {
      k1 = Normalize(r)*adt(Length(r)); //First step.
      r2 = r + oldvel*hh;
      k2 = Normalize(r2)*adt(Length(r2));  //Second step.
      r2 = r + oldveldt + k1*dt;
      k3 = Normalize(r2)*adt(Length(r2));

      v += (k1 + k2*4 + k3)*(1/6.0);
      x += (k1 + k2*2)*(1/6.0)*dt;
    }
  }
  p.vel[i] += v;
  p.pos[i] += oldveldt + x;
}


void ComputeRK3()
{
  int i;

  for (i = 1; i <= nsteps; i++)
    if (!p.finished[0])
    {
      rk3(0);
      if (p.pos[0].z > zmax)
      {
        p.finished[0] = true;
        stoped = true;
      }
    }
  if (!p.finished[0])
  {
    Track[jj] = p.pos[0];
    jj += 1;
  }
}




void rk4(int i)
{
  int j;
  VECTOR r, r2;
  VECTOR k1, k2, k3, k4;

  VECTOR v,x;
  VECTOR oldpos, oldvel, oldvelhh, oldveldt;
  REAL modul_r;


  Null(x);
  Null(v);

  oldpos = p.pos[i];
  oldvel = p.vel[i];
  oldvelhh = oldvel*hh;
  oldveldt = oldvel*dt;

  for (j = 0; j <= 349; j++)
  {
    r = oldpos - A.pos[j];
    modul_r = Length(r);
    if (modul_r < radius)
    {
      k1 = Normalize(r)*adt(Length(r)); //First step.
      r2 = r + oldvelhh;
      k2 = Normalize(r2)*adt(Length(r2));  //Second step.
      r2 = r + oldvelhh + k1*0.5*hh;
      k3 = Normalize(r2)*adt(Length(r2));
      r2 = r + oldveldt + k2*hh;
      k4 = Normalize(r2)*adt(Length(r2));

      v += (k1 + k2*2 + k3*2 + k4)*(1/6.0);
      x += (k1 + k2 + k3)*h6;
    }
  }
  p.vel[i] += v;
  p.pos[i] += oldveldt + x;
}


void ComputeRK4()
{
  int i;

  for (i = 1; i <= nsteps; i++)
  {
    if (!p.finished[0])
    {
      rk4(0);
      if (p.pos[0].z > zmax)
      {
        p.finished[0] = true;
        stoped = true;
      }
    }
  }
  if (!p.finished[0])
  {
    Track[jj] = p.pos[0];
    jj += 1;
  }

}



void rkck(int i, VECTOR *v, VECTOR *x)
{
  int j, ii, jj;
  VECTOR r, r2;

  VECTOR err1, err2;
  VECTOR oldpos, oldvel, oldveldt;
  REAL modul_r;

  VECTOR kv[6], kr[6], krinc[6];

  Null(*x);
  Null(*v);
  Null(err1);
  Null(err2);

  oldpos = p.pos[i];
  oldvel = p.vel[i];
  oldveldt = oldvel*h;

  kr[0] = oldveldt;
  Null(krinc[0]);

  for (j = 0; j <= 349; j++)
  {
     r = oldpos - A.pos[j];
     modul_r = Length(r);
     if (modul_r < radius)
     {
       kv[0] = Normalize(r)*ah(Length(r));
       for(ii = 1; ii <= 5; ii++)
       {
         r2 = r;
         for(jj = 0; jj <= ii-1; jj++) r2 += kr[jj]*b[ii][jj];
         kv[ii] = Normalize(r2)*ah(Length(r2));
         Null(krinc[ii]);
         for(jj = 0; jj <= ii-1; jj++) krinc[ii] += kv[jj]*b[ii][jj];
         krinc[ii] = krinc[ii]*h;
         kr[ii] = kr[0] + krinc[ii];
       }

       *v += kv[0]*c1 + kv[2]*c3 + kv[3]*c4 + kv[5]*c6;
       *x += krinc[2]*c3 + krinc[3]*c4 + krinc[5]*c6;

       err1 += kv[0]*dc1 + kv[2]*dc3 + kv[3]*dc4 + kv[5]*dc6;
       err2 += krinc[2]*dc3 + krinc[3]*dc4 + krinc[5]*dc6;
     }
  }
  *x += kr[0];


  e1 = Length(err1);
  e2 = Length(err2);
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4


//step-monitoring proc

void rksm (int i)
{


  REAL htemp, hnext, errmax, eps, yscal, tnew;
  VECTOR v, x;


//  eps = 3000*4.0e-8;
//  yscal = h*Length(p.vel[0]);
  eps = 4.0e-7;
  yscal = 4.0;


 while( t-tsav < Tsav ){
  {
    for(;;)
    {

      rkck(i, &v, &x);
      errmax = e2/yscal;
      errmax /= eps;
      if (errmax <= 1.0) {
        nok++;
        p.vel[i] += v;
        p.pos[i] += x;
        break;
      }
      nbad++;
      htemp = SAFETY*h*pow(errmax, PSHRNK);
      h = FMAX(htemp, 0.1*h);
      tnew = t + h;
      if (tnew == t) AfxMessageBox("Stepsize underflow in rksm!");
    }
    if (errmax > ERRCON) hnext = SAFETY*h*pow(errmax, PGROW);
    else hnext = 5.0*h;
    t += h;
  }
  h = hnext;
 }
 tsav = t;
 if (!p.finished[0])
 {
  Track[jj] = p.pos[0];
  jj += 1;
 }
}




void ComputeRKCK()
{
//  int i;

  if (!p.finished[0])
  {
    rksm(0);
    if (p.pos[0].z > zmax)
    {
      p.finished[0] = true;
      stoped = true;
    }
  }
}



//===========================================================
//========== Вспомогательный класс
class CPoint3D
{
public:
	//====== Координаты точки
	float x;
	float y;
	float z;

	//====== Набор конструкторов
	CPoint3D () 
	{
 	x = y = z = 0.f;
 }

	CPoint3D (float c1, float c2, float c3)
	{
		x = c1;
		z = c2;
		y = c3;
	}
	//====== Операция присвоения
	CPoint3D& operator= (const CPoint3D& pt)
	{
		x = pt.x;	z = pt.z;	y = pt.y;
		return *this;
	}

	//====== Конструктор копирования
	CPoint3D (const CPoint3D& pt)
	{
		*this = pt;
	}
};



//===================================================
//===================================================



	CListBox* pBox;
	






class CObjectWnd : public CDialog
{

  public:


	CObjectWnd();	
	CListCtrl* pProp;
	afx_msg BOOL OnSetCursor( CWnd* pWnd, UINT nHitTest, UINT message );
	afx_msg void OnCancel();
	//afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	//afx_msg void OnPaint();
	afx_msg void OnClick();
	virtual  BOOL OnInitDialog( );
	void ShowEn();

  protected:	
    DECLARE_MESSAGE_MAP();


	
};

BEGIN_MESSAGE_MAP(CObjectWnd,CDialog)
 ON_WM_SETCURSOR()
 ON_WM_CLOSE()
// ON_WM_SYSCOMMAND()
// ON_WM_PAINT()
 ON_LBN_SELCHANGE (IDC_LIST1, OnClick)
END_MESSAGE_MAP()


CObjectWnd * ObjectWnd;





void CObjectWnd::ShowEn()
{

  REAL E, AbsErrorE, RelErrorE, AvRelErrorE;

  E = Energy();
  AbsErrorE = E0-E;
  RelErrorE = fabs(AbsErrorE/E0);
  SummREE += RelErrorE;
  AvRelErrorE = SummREE/(jj-2);




	pProp = (CListCtrl*) GetDlgItem(IDC_LIST2);
  	char s[200];
	sprintf(s, "%e", E0);
	pProp->SetItemText(0,1,s);
	sprintf(s, "%e", Ek);
	pProp->SetItemText(1,1,s);
	sprintf(s, "%e", Ep);
	pProp->SetItemText(2,1,s);









/*  Form_EnergyList->Label2->Caption = AnsiString(double(Ek));
  Form_EnergyList->Label4->Caption = AnsiString(double(Ep));
  Form_EnergyList->Label6->Caption = AnsiString(double(E0));

  Form_EnergyList->Label8->Caption = AnsiString(double(E));
  Form_EnergyList->Label10->Caption = AnsiString(double(AbsErrorE));
  Form_EnergyList->Label12->Caption = AnsiString(double(RelErrorE));
  Form_EnergyList->Label14->Caption = AnsiString(double(AvRelErrorE));

  Form_EnergyList->Label18->Caption = AnsiString(double(e1));
  Form_EnergyList->Label16->Caption = AnsiString(double(e2));

  Form_EnergyList->Label20->Caption = AnsiString(double(dt));
  Form_EnergyList->Label22->Caption = AnsiString(double(h));


  Form_EnergyList->Label27->Caption = AnsiString(double(t));


  Form_EnergyList->Label24->Caption = AnsiString(nok);
  Form_EnergyList->Label26->Caption = AnsiString(nbad);*/



}



void CObjectWnd::OnClick()
{
	
	pBox = (CListBox*) GetDlgItem(IDC_LIST1);
	int index = pBox->GetCurSel();

	//char s[200];
	//sprintf(s, "%i", index);


	pProp = (CListCtrl*) GetDlgItem(IDC_LIST2);
	pProp->SetItemText(index,1,"111");




}



BOOL CObjectWnd::OnSetCursor( CWnd* pWnd, UINT nHitTest, UINT message )
{


	this->SetFocus();
	return CDialog::OnSetCursor(pWnd, nHitTest, message);


}




BOOL CObjectWnd::OnInitDialog()
{

  CDialog::OnInitDialog();

	HFONT hFont;
	hFont = CreateFont(
		14,                     // nHeight
		6,                      // nWidth
		0,                      // nEscapement
		0,                      // nOrientation
		FW_NORMAL,                // nWeight
		0,                      // bItalic
		0,                      // bUnderline
		0,                      // cStrikeOut
		GREEK_CHARSET,          // nCharSet
		OUT_DEFAULT_PRECIS,     // nOutPrecision
		CLIP_DEFAULT_PRECIS,    // nClipPrecision
		PROOF_QUALITY,          // nQuality
		34,						// nPitchAndFamily
		"Arial");



  	pBox = (CListBox*) GetDlgItem(IDC_LIST1);
	pBox->SendMessage(WM_SETFONT, (WPARAM)hFont, TRUE);	
	pBox->AddString("p[0]");
	pBox->AddString("p[1]");
	pBox->AddString("p[2]");




	


	pProp = (CListCtrl*) GetDlgItem(IDC_LIST2);
	pProp->SendMessage(WM_SETFONT, (WPARAM)hFont, TRUE);
	pProp->InsertColumn(0, "Property", LVCFMT_LEFT, 60, -1);
	pProp->InsertColumn(1, "Value", LVCFMT_LEFT, 124, -1);
	pProp->InsertItem(0, 0);
	pProp->SetItemText(0,0, "Eo");
	pProp->InsertItem(1, 0);
	pProp->SetItemText(1,0, "Ek");
	pProp->InsertItem(2, 0);
	pProp->SetItemText(2,0, "Ep");
	pProp->InsertItem(3, 0);
	pProp->SetItemText(3,0, "Ek+Ep");
	pProp->InsertItem(4, 0);
	pProp->SetItemText(4,0, "ДE");
	pProp->InsertItem(5, 0);
	pProp->SetItemText(5,0, "дE");
	



  return TRUE; 

}

CObjectWnd::CObjectWnd()
{


	Create(IDD_DIALOG4, this);


}





void CObjectWnd::OnCancel()
{  

   DestroyWindow();

}






//====================================================
//====================================================


class CViewWnd : public CDialog
{

  public:

	CViewWnd ();	
	afx_msg BOOL OnSetCursor ( CWnd* pWnd, UINT nHitTest, UINT message );
	afx_msg void OnSysCommand ( UINT nID, LPARAM lParam );

	afx_msg void OnPaint ();
	afx_msg void OnCancel ();
	afx_msg void OnSize (UINT nType, int cx, int cy);

	afx_msg void OnLButtonDown( UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp( UINT nFlags, CPoint point);
	afx_msg void OnMouseMove( UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel( UINT nFlags, short zDelta, CPoint pt );

	void OnXY ();
	void OnYZ ();
	void OnZX ();
	void OnOrtho ();
	void OnPerspective ();

	void CreateGL();
	void ConfigGL(void);
	void CloseGL();

	GLfloat	m_AngleX;
	GLfloat	m_AngleY;

	GLfloat	m_ax;
	GLfloat	m_ay;
	GLfloat	m_az;

	GLfloat	m_xTrans;
	GLfloat	m_yTrans;
	GLfloat	m_zTrans;

	CPoint		m_pt;				// Текущая позиция мыши
	bool		m_bCaptured;		// Признак захвата мыши

	GLfloat	m_dx;				// Квант смещения вдоль X
	GLfloat	m_dy;				// Квант смещения вдоль Y

	UINT		m_cx;
	UINT		m_cy;

	BOOL m_Perspective;

	HDC hDC;
	HGLRC hRC;

	
	

  protected:	
    DECLARE_MESSAGE_MAP();
	virtual void OnOK();
	
};

BEGIN_MESSAGE_MAP(CViewWnd,CDialog)
 ON_WM_SETCURSOR( )
 ON_WM_SYSCOMMAND( )
 ON_WM_PAINT()
 ON_WM_CLOSE()
 ON_WM_SIZE()
 ON_WM_LBUTTONDOWN()
 ON_WM_LBUTTONUP()
 ON_WM_MOUSEMOVE()
 ON_WM_MOUSEWHEEL()
END_MESSAGE_MAP()





CViewWnd*	p_ViewWnd;
CWnd*	p_MainWnd;






void CViewWnd::OnOK() { }


CViewWnd::CViewWnd()
{
	

	Create(IDD_DIALOG3, p_MainWnd);
	CreateGL();
	ConfigGL();


}

void CViewWnd::CreateGL()
{
  
   
  hDC = ::GetDC(GetSafeHwnd());
  
	PIXELFORMATDESCRIPTOR pfd =	
	{
		sizeof(PIXELFORMATDESCRIPTOR),
		1,							
		PFD_DRAW_TO_WINDOW |	
		PFD_SUPPORT_OPENGL |	
		PFD_DOUBLEBUFFER,		
		PFD_TYPE_RGBA,			
		24, 					
		 						
		24,	0,					
		24,	0,					
		24,	0,					
		24,	0,					
		0,						
								
		0,			 			
		0,			 			
		0,						
		0,						
		32, 					
		0,						
		0,						
		0,				 		
		0,						
		0,						
		0,						
		0
	};

	int i = ChoosePixelFormat( hDC, &pfd );
	SetPixelFormat( hDC, i, &pfd );
	hRC = wglCreateContext( hDC );
	wglMakeCurrent( hDC, hRC );

  
}


void CViewWnd::ConfigGL(void)
{


	glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
	glShadeModel( GL_SMOOTH );
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glEnable( GL_LINE_SMOOTH );


	m_AngleX = 0.f;
	m_AngleY = 0.f;

	m_xTrans = -10.f;
	m_yTrans = 0.f;
	m_zTrans = 0.f;

	m_ax = 0;
	m_ay = 90;
	m_az = 90;

	m_Perspective = FALSE;

	glEnable (GL_POINT_SMOOTH);
	glHint (GL_POINT_SMOOTH_HINT, GL_DONT_CARE);	


	glEnable(GL_DEPTH_TEST);
	glEnable(GL_ALPHA_TEST);

}



void CViewWnd::OnSize(UINT nType, int cx, int cy)
{


	CDialog::OnSize(nType, cx, cy);

	int inc1, inc2, ind;
	int i;

	m_cx = cx;
	m_cy = cy;

	glViewport(0, 0, cx, cy);
	glPointSize(float(cx/200.0));
	

	glNewList(1, GL_COMPILE);

		glBegin(GL_POINTS);
		for(i = 0; i <= 174; i++)
		{

			glColor4d(0.9,0.5,0.0,0.9);
			glVertex3d(A.pos[i].x, A.pos[i].y, A.pos[i].z);
			glColor4d(0.5,0.5,1.0,0.9);
			glVertex3d(A.pos[i+175].x, A.pos[i+175].y, A.pos[i+175].z);

		}
		glEnd();


		//glLineWidth(float(cx/333.0));
		glLineWidth( 1.5 );
		glBegin( GL_LINES );
		for( i = 0; i <= 174; i++ )
		{

			inc1 = 0, inc2 = 0;
			glColor4d( 0.9,0.9,0.9, 0.3 );
			glVertex3d( A.pos[i].x, A.pos[i].y, A.pos[i].z );
			ind = 175 + (i+7)%175;
			glVertex3d( A.pos[ind].x, A.pos[ind].y, A.pos[ind].z );

			if ( ((i/7-2)%5==0) & ((i+1)%7==0) )
			{

				glVertex3d( A.pos[i].x, A.pos[i].y, A.pos[i].z );
				ind = i+175;
				glVertex3d( A.pos[ind].x, A.pos[ind].y, A.pos[ind].z );

			} else {

				if( !(( ((i/7)%5==0) | ((i/7-3)%5==0) ) & ((i)%7==0) ) )
				{

					if ((i/7-2)%5 == 0) inc2=1;
					if ( ((i/7)%5 == 0) | ((i/7-3)%5 == 0) ) inc1=-1;

					ind = i + 175 + inc1;
					glVertex3d( A.pos[i].x, A.pos[i].y, A.pos[i].z );
					glVertex3d( A.pos[ind].x, A.pos[ind].y, A.pos[ind].z );

				}
				ind = 175 + (175 + i - 7 + inc2)%175;
				glVertex3d( A.pos[i].x, A.pos[i].y, A.pos[i].z );
				glVertex3d( A.pos[ind].x, A.pos[ind].y, A.pos[ind].z );

			}
		}
		glEnd( );


		glBegin( GL_LINES );
			glColor4d( 0.9, 0.4, 0, 0.5 );
		    glVertex3d( -5.0, -5.0, 0.0 );
			glVertex3d( 5.0, -5.0, 0.0 );

			glColor4d( 0.0, 0.9, 0.5, 0.5 );
			glVertex3d( -5.0, -5.0, 0.0 );
			glVertex3d( -5.0, 5.0, 0.0 );

			glColor4d( 0.0, 0.5, 0.9, 0.5 );
		    glVertex3d( -5.0, -5.0, 0.0 );
			glVertex3d( -5.0, -5.0, 20.0 );
		glEnd();
	glEndList();

  
	Invalidate(FALSE);



}


void CViewWnd::CloseGL()
{

	wglMakeCurrent( NULL, NULL );
	wglDeleteContext( hRC );
	hRC = 0;

}







void CViewWnd::OnPaint()
{

	

	double dAspect = m_cx<=m_cy ? double(m_cy)/m_cx : double(m_cx)/m_cy;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();

	if (m_Perspective == TRUE)
	{

		gluPerspective(23.0, dAspect, 0.01, 10000);
		gluLookAt(0, 0, 5, 0, 0, 0,  0,  1, 0);
		glScaled( scale, scale, scale );

	} else
		glScaled( 0.7*scale, scale, scale );



	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glTranslatef(m_xTrans, m_yTrans, m_zTrans);
	glRotatef (m_AngleX+m_ax, 1.0f, 0.0f, 0.0f );
	glRotatef (m_AngleY+m_ay, 0.0f, 1.0f, 0.0f);
	glRotatef (m_az, 0.0f, 0.0f, 1.0f);


	//glLineWidth(2);
	glBegin( GL_LINES );
	glColor3d(0.9,0.9,0.9);
	int i;
	for (i = 0; i < int(jj-1); i++)
	{
		glVertex3d(double(Track[i].x), double(Track[i].y), double(Track[i].z));
		glVertex3d(double(Track[i+1].x), double(Track[i+1].y), double(Track[i+1].z));
	}
	glEnd();
	//glLineWidth(1.0);


	glBegin(GL_POINTS);
		glColor3d(0.9,0.9,0.1);
		glVertex3d(p.pos[0].x, p.pos[0].y, p.pos[0].z);
	glEnd ();	
	
	glCallList(1);

	SwapBuffers(hDC);

	CDialog::OnPaint();


}


void CViewWnd::OnLButtonDown( UINT nFlags, CPoint point)
{

	SetCapture();
	m_bCaptured = true;
	m_pt = point;
	

}

void CViewWnd::OnLButtonUp( UINT nFlags, CPoint point)
{

	if (m_bCaptured == true)
	{

		m_bCaptured = false;
		ReleaseCapture();

	}

}


BOOL CViewWnd::OnSetCursor( CWnd* pWnd, UINT nHitTest, UINT message )
{

	this->SetFocus();
	return CDialog::OnSetCursor(pWnd, nHitTest, message);

}

void CViewWnd::OnMouseMove( UINT nFlags, CPoint point)
{

	
	if (m_bCaptured == true)
	{

		m_dy = float(point.y - m_pt.y)/5.f;
		m_dx = float(point.x - m_pt.x)/5.f;

		m_AngleX += m_dy;
		m_AngleY += m_dx;
		

		m_pt = point;

		Invalidate(FALSE);

	}

}

void CViewWnd::OnCancel()
{  

	timer = TimerNonInteractive;
	this->ShowWindow(SW_MINIMIZE);

}

void CViewWnd::OnSysCommand( UINT nID, LPARAM lParam )
{

	if ( nID == SC_RESTORE )
	{
		timer = TimerInteractive;
	}

	CDialog::OnSysCommand(nID, lParam);

}


void CViewWnd::OnXY ()
{

	m_ax = 0;
	m_ay = 180;
	m_az = 90;

	m_xTrans = 0.f;
	m_yTrans = 0.f;
	m_zTrans = 0.f;

	m_AngleX = 0;
	m_AngleY = 0;

	Invalidate(FALSE);

}

void CViewWnd::OnYZ ()
{

	m_ax = 0;
	m_ay = 90;
	m_az = 180;

	m_xTrans = -10.f;
	m_yTrans = 0.f;
	m_zTrans = 0.f;

	m_AngleX = 0;
	m_AngleY = 0;

	Invalidate(FALSE);
	
}

void CViewWnd::OnZX ()
{

	m_ax = 0;
	m_ay = 90;
	m_az = 90;

	m_xTrans = -10.f;
	m_yTrans = 0.f;
	m_zTrans = 0.f;

	m_AngleX = 0;
	m_AngleY = 0;

	Invalidate(FALSE);

}

void CViewWnd::OnOrtho ()
{

	m_Perspective = FALSE;
	Invalidate(FALSE);

}

void CViewWnd::OnPerspective ()
{

	m_Perspective = TRUE;
	Invalidate(FALSE);

}

BOOL CViewWnd::OnMouseWheel( UINT nFlags, short zDelta, CPoint pt )
{

	scale += scale/100*zDelta/120;
	Invalidate(FALSE);
	return TRUE;

}














//=================================================
//=================================================



class CAbout : public CDialog
{

	public:
	CAbout(CWnd* pParent = NULL);   // standard constructor
	enum { IDD = IDD_DIALOG2 };
	protected:	
	virtual void OnOK();

};



void CAbout::OnOK()
{

	CDialog::OnOK();
	p_ViewWnd->EnableWindow(TRUE);
	ObjectWnd->EnableWindow(TRUE);

}

CAbout::CAbout(CWnd* pParent)
	:CDialog(CAbout::IDD, pParent)
{

	p_ViewWnd->EnableWindow(FALSE);
	ObjectWnd->EnableWindow(FALSE);

	
}





class CMainWnd : public CDialog
{

  public:



	CMainWnd();	

	afx_msg BOOL OnSetCursor( CWnd* pWnd, UINT nHitTest, UINT message );
	afx_msg void OnGetMinMaxInfo(MINMAXINFO* lpMMI);
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg void OnCancel();
	afx_msg void OnAbout();
	afx_msg void OnExit();
	afx_msg void OnPerspective();
	afx_msg void OnOrtho();
	afx_msg void OnXY();
	afx_msg void OnYZ();
	afx_msg void OnZX();
	virtual void DoDataExchange(CDataExchange* pDX);
  protected:	
    DECLARE_MESSAGE_MAP();


	
};

BEGIN_MESSAGE_MAP(CMainWnd,CDialog)
 ON_WM_SETCURSOR()
 ON_WM_CLOSE()
 ON_BN_CLICKED(ID_ABOUT, OnAbout)
 ON_BN_CLICKED(ID_EXIT, OnExit)
 ON_BN_CLICKED(ID_PERSPECTVE, OnPerspective)
 ON_BN_CLICKED(ID_ORTHO, OnOrtho)
 ON_BN_CLICKED(ID_XY, OnXY)
 ON_WM_GETMINMAXINFO()
 ON_WM_SYSCOMMAND()
 ON_WM_PAINT()
 ON_BN_CLICKED(ID_YZ, OnYZ)
 ON_BN_CLICKED(ID_ZX, OnZX)
END_MESSAGE_MAP()





BOOL CMainWnd::OnSetCursor( CWnd* pWnd, UINT nHitTest, UINT message )
{


	this->SetFocus();
	return CDialog::OnSetCursor(pWnd, nHitTest, message);


}





CMainWnd::CMainWnd()
{

	Create(IDD_DIALOG1, this);

	menu.LoadMenu(IDR_MENU1);
	this->SetMenu(&menu);

	Bar.Create(WS_VISIBLE|WS_CHILD|TBSTYLE_FLAT, CRect(0,0,::GetSystemMetrics(SM_CXSCREEN),0), this,1);
	TBBUTTON b[6];
	b[0].iString = Bar.AddStrings(_T("Perspective"));
	b[0].fsState = TBSTATE_ENABLED;
	b[0].fsStyle = TBSTYLE_GROUP;
	b[0].iBitmap = -1;
	b[0].idCommand = ID_PERSPECTVE;

	b[1].iString = Bar.AddStrings(_T("Ortho"));
	b[1].fsState = TBSTATE_ENABLED;
	b[1].fsStyle = TBSTYLE_GROUP;
	b[1].iBitmap = -1;
	b[1].idCommand = ID_ORTHO;

	b[5].fsState = TBSTATE_ENABLED;
	b[5].fsStyle = TBSTYLE_SEP;
	b[5].iBitmap = -1;

	b[2].iString = Bar.AddStrings(_T("XY"));
	b[2].fsState = TBSTATE_ENABLED;
	b[2].fsStyle = TBSTYLE_GROUP;
	b[2].iBitmap = -1;
	b[2].idCommand = ID_XY;

	b[3].iString = Bar.AddStrings(_T("YZ"));
	b[3].fsState = TBSTATE_ENABLED;
	b[3].fsStyle = TBSTYLE_GROUP;
	b[3].iBitmap = -1;
	b[3].idCommand = ID_YZ;

	b[4].iString = Bar.AddStrings(_T("ZX"));
	b[4].fsState = TBSTATE_ENABLED;
	b[4].fsStyle = TBSTYLE_GROUP;
	b[4].iBitmap = -1;
	b[4].idCommand = ID_ZX;


	Bar.InsertButton(0, &b[0]);
	Bar.InsertButton(1, &b[1]);
	Bar.InsertButton(2, &b[5]);
	Bar.InsertButton(3, &b[2]);
	Bar.InsertButton(4, &b[3]);
	Bar.InsertButton(5, &b[4]);

	Bar.SetButtonSize(CSize(10,40));	
	Bar.SetBitmapSize(CSize(0,0));
	Bar.SetButtonWidth(0,200);
	Bar.AutoSize();

	



}


void CMainWnd::OnGetMinMaxInfo(MINMAXINFO* lpMMI)
{

	lpMMI->ptMaxSize.y = 67;

}

void CMainWnd::OnSysCommand(UINT nID, LPARAM lParam)
{
	
	CWnd::OnSysCommand(nID, lParam);

}


void CMainWnd::OnPaint()
{
    
	CDialog::OnPaint();

}

void CMainWnd::OnCancel()
{  

   DestroyWindow();
   


}

void CMainWnd::OnAbout()
{  

	CAbout aboutDlg;
	aboutDlg.DoModal();	

}

void CMainWnd::OnExit()
{  

	DestroyWindow();

}

void CMainWnd::DoDataExchange(CDataExchange* pDX)
{

	CDialog::DoDataExchange(pDX);

}

void CMainWnd::OnPerspective()
{

	p_ViewWnd->OnPerspective();
	
}

void CMainWnd::OnOrtho()
{  

	p_ViewWnd->OnOrtho();
}

void CMainWnd::OnXY()
{  

	p_ViewWnd->OnXY();
	
}

void CMainWnd::OnYZ()
{  

	p_ViewWnd->OnYZ();

}

void CMainWnd::OnZX()
{  

	p_ViewWnd->OnZX();

}




//====================================================
//====================================================




class CMyApp : public CWinApp
{

  public:

	CMyApp();			
	CViewWnd *m_pViewWnd;
	int LoadData();
	virtual BOOL InitInstance();
	virtual BOOL PreTranslateMessage ( MSG* pMsg );
	int CMyApp::ExitInstance();
	
	
  protected:	
    DECLARE_MESSAGE_MAP();
	long double ReadFloat(LPCTSTR section, LPCTSTR entry);
	int ReadInt(LPCTSTR section, LPCTSTR entry);
	CString ReadString(LPCTSTR section, LPCTSTR entry);

};

BEGIN_MESSAGE_MAP(CMyApp, CWinApp)
END_MESSAGE_MAP()



CMyApp theApp;






CMyApp::CMyApp()
{

}

long double CMyApp::ReadFloat(LPCTSTR section, LPCTSTR entry)
{

	CString STemp = "";
	long double LD = 0.0;
	STemp = AfxGetApp()->GetProfileString(section, entry, 0);
	sscanf((const char *)STemp, "%Lf", &LD);
	return LD;

	

}

int CMyApp::ReadInt(LPCTSTR section, LPCTSTR entry)
{

	CString STemp = "";
	int Int = 0;
	STemp = AfxGetApp()->GetProfileString(section, entry, 0);
	sscanf((const char *)STemp, "%i", &Int);
	return Int;

}

CString CMyApp::ReadString(LPCTSTR section, LPCTSTR entry)
{

	CString STemp = "";
	STemp = AfxGetApp()->GetProfileString(section, entry, 0);
	return STemp;

}


void Parse(char fname[255], DWORD nulloffset, int step)
{

  FILE *f;
  int unsigned i, j;
  char s[300] = "";
  char temp[300] = "";
  REAL number;
  DWORD adress;
  REAL * r_pointer;

  f = fopen(fname, "r");
  adress = nulloffset;
  while (!feof(f))
  {
    strset(s, 0);
    fscanf(f, "%s", s);
    if (strlen(s)>4)
    {
      strncpy(s+strlen(s), ",", 2);
      j = 0;
      for (i=1; i<=strlen(s); i++)
      {
        if (s[i] == ',')
        {
          strset(temp, 0);
          strncpy(temp, s+j, i-j);
          sscanf(temp, "%Lf", &number);
		  
		  //r_pointer = adress;
		  _asm
		  {
			  mov eax, adress
			  mov r_pointer, eax
		  }
          *r_pointer = number;


          adress += step;
          j = i + 1;
        }
      }
    }
  }
  fclose(f);
}


int CMyApp::LoadData()
{


	DWORD step = offsetof(ATOMS, pos[1].x)-offsetof(ATOMS, pos[0].x);
	DWORD nulloffset = DWORD(&A)+offsetof(ATOMS, pos[0].x);
	Parse("x1.csv", nulloffset, step);
	nulloffset = DWORD(&A)+offsetof(ATOMS, pos[0].y);
	Parse("y1.csv", nulloffset, step);
	nulloffset = DWORD(&A)+offsetof(ATOMS, pos[0].z);
	Parse("z1.csv", nulloffset, step);
	nulloffset = DWORD(&A)+offsetof(ATOMS, pos[175].x);
	Parse("x2.csv", nulloffset, step);
	nulloffset = DWORD(&A)+offsetof(ATOMS, pos[175].y);
	Parse("y2.csv", nulloffset, step);
	nulloffset = DWORD(&A)+offsetof(ATOMS, pos[175].z);
	Parse("z2.csv", nulloffset, step);


	
	char c_dirname[255];
	GetCurrentDirectory(255, c_dirname);
	this->m_pszProfileName = _tcscat(c_dirname, "\\particle.ini");	
	
	scale =		ReadFloat("ENGINE", "scale");
	Tx =		ReadFloat("ENGINE", "T.x");
	Ty =		ReadFloat("ENGINE", "T.y");
	Tz =		ReadFloat("ENGINE", "T.z");
	dt =		ReadFloat("ENGINE", "dt");
	Tsav =		ReadFloat("ENGINE", "Tsav");
	radius =	ReadFloat("ENGINE", "radius");
	zmax =		ReadFloat("ENGINE", "zmax");
	nsteps =	ReadInt("ENGINE", "nsteps");



	if (ReadString("ENGINE", "interactive") == "yes")
	{
		timer = TimerInteractive;
	}
	else
	{
		timer = TimerNonInteractive;
	}


	CString method = ReadString("ENGINE", "method");
	Compute = ComputeEuler;
	char newtitle[100];
	CString title;	
	p_MainWnd->GetWindowText(title);
	sprintf(newtitle, "%s using Euler method", title);
	
	if (method=="RK2")
	{
	    initRK();
		Compute = ComputeRK2;
		sprintf(newtitle, "%s using RK2 method", title);
	}
	if (method=="RK3")
	{
		initRK();
		Compute = ComputeRK3;
		sprintf(newtitle, "%s using RK3 method", title);
	}
	if (method=="RK4")
	{
		initRK();
		Compute = ComputeRK4;
		sprintf(newtitle, "%s using RK4 method", title);
	}
	if (method=="RKCK")
	{
	    initRK();
		Compute = ComputeRKCK;
		sprintf(newtitle, "%s using RK-Cash-Karp method", title);
	}

	p_MainWnd->SetWindowText(newtitle);


	p.charge =   (float)ReadInt("PARTICLE", "q");
	p.mass =     ReadFloat("PARTICLE", "m");
	p.vel[0].x = ReadFloat("PARTICLE", "vel.x")*iatf;
	p.vel[0].y = ReadFloat("PARTICLE", "vel.y")*iatf;
	p.vel[0].z = ReadFloat("PARTICLE", "vel.z")*iatf;

	p.pos[0].x = ReadFloat("PARTICLE", "pos.x");
	p.pos[0].y = ReadFloat("PARTICLE", "pos.y");
	p.pos[0].z = ReadFloat("PARTICLE", "pos.z");

	p.charge = -p.charge;
	_C = 80750.750385769996443*p.charge/p.mass;
	C = _C*dt;
	Cek = halfatf2*p.mass;
	Cpot = -5.3637847231117074279e-17*p.charge;




	//char CTemp[255];
	//sprint(CTemp, "%12.12e", dt);
	//sprintf(CTemp, "%i", nsteps);
	//AfxMessageBox(CTemp,0,0);


	return 0;
}









void TimerInteractive()
{


	
	ObjectWnd->ShowEn();
	
	p_ViewWnd->Invalidate(FALSE);


}

void TimerNonInteractive()
{


	
	ObjectWnd->ShowEn();

	
}



BOOL CMyApp::PreTranslateMessage ( MSG* pMsg )
{

	if( pMsg->message == WM_NEWDATA )
    {

		timer();
		ready = TRUE;
        return TRUE;

    }
  
    
	
	BOOL rvalue = CWinApp::PreTranslateMessage ( pMsg );
	ready = TRUE;
	return rvalue;
	

}





UINT ThreadCode(LPVOID pParam)
{

	while(!stoped)
	{
		Sleep(1);
		if (ready == TRUE)
		{
			Compute();
			ready = FALSE;
			PostThreadMessage(thread_id, WM_NEWDATA, 0, 0);
		}
	}
	

	return 0;

}


BOOL CMyApp::InitInstance() 
{


	Enable3dControls();

	p_MainWnd = m_pMainWnd = new CMainWnd();
	m_pMainWnd->MoveWindow(0, 0, ::GetSystemMetrics(SM_CXSCREEN), 67);
	m_pMainWnd->ShowWindow(SW_SHOW);
	m_pMainWnd->UpdateWindow();


	LoadData();
	
		
	
	p_ViewWnd = m_pViewWnd = new CViewWnd();
	m_pViewWnd->MoveWindow(1, 67, 576, 442);
	m_pViewWnd->ShowWindow(SW_SHOW);
	m_pViewWnd->UpdateWindow();
	m_pMainWnd->SetFocus();

	ObjectWnd = new CObjectWnd();
	ObjectWnd->MoveWindow(577, 67, 200, 442);
	ObjectWnd->ShowWindow(SW_SHOW);
	ObjectWnd->UpdateWindow();
	ObjectWnd->SetFocus();
	
	
	


	thread_id = this->m_nThreadID;
	CWinThread *computer;	
	computer = new CWinThread();
	//computer->SetThreadPriority(THREAD_PRIORITY_NORMAL);	
	ready = TRUE;
	AfxBeginThread(ThreadCode, computer);


	
	

	return TRUE;

};


int CMyApp::ExitInstance() 
{

	return CWinApp::ExitInstance();

}
