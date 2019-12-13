#include <cppunit/extensions/HelperMacros.h>

#include <MathTools.hh>

class Tools_ScaleGraphAndMove : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE( Tools_ScaleGraphAndMove );
	CPPUNIT_TEST( ScaleTest );
	CPPUNIT_TEST( MoveTest );
    CPPUNIT_TEST( ScaleTGraphErrorsTest );
    CPPUNIT_TEST( MoveTGraphErrorsTest );
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
	void ScaleTest();
	void MoveTest();
    void ScaleTGraphErrorsTest();
    void MoveTGraphErrorsTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION( Tools_ScaleGraphAndMove );

void Tools_ScaleGraphAndMove::setUp() {
}

void Tools_ScaleGraphAndMove::tearDown() {
}

void Tools_ScaleGraphAndMove::ScaleTest() {
	const size_t n = 5;
	double x[n] = { 1, 2, 3, 4, 5 };
	double y[n] = { 5, 4, 3, 2, 1 };

	TGraph *gr = new TGraph(n, x, y);

	MathTools::ScaleGraphAndMove(gr, 2.0);

	for (uint i = 0; i < n; ++i) {
		double _x, _y;
		gr->GetPoint(i, _x, _y);
		CPPUNIT_ASSERT_DOUBLES_EQUAL( y[i] * 2.0, _y, 0.00001);
	}
	delete gr;
}

void Tools_ScaleGraphAndMove::MoveTest() {
	const size_t n = 5;
	double x[n] = { 1, 2, 3, 4, 5 };
	double y[n] = { 5, 4, 3, 2, 1 };

	TGraph *gr = new TGraph(n, x, y);

	MathTools::ScaleGraphAndMove(gr, 1.0, -5.0);

	for (uint i = 0; i < n; ++i) {
		double _x, _y;
		gr->GetPoint(i, _x, _y);
		CPPUNIT_ASSERT_DOUBLES_EQUAL( x[i] - 5.0, _x, 0.00001);
	}
	delete gr;
}

void Tools_ScaleGraphAndMove::ScaleTGraphErrorsTest() {
    const size_t n = 5;
    double x[n]  = { 1, 2, 3, 4, 5 };
    double ex[n] = { 0.1, 0.1, 0.1, 0.1, 0.1 };
	double y[n]  = { 5, 4, 3, 2, 1 };
    double ey[n] = { 0.1, 0.1, 0.1, 0.1, 0.1 };
    
    TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);
    
    MathTools::ScaleGraphAndMove(gr, 2.0);
    
    for(uint i = 0; i<n; ++i){
        double _x, _y;
        double _ey = gr->GetErrorY(i);
        gr->GetPoint(i, _x, _y);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( y[i] * 2.0, _y, 0.00001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( ey[i] * 2.0, _ey, 0.00001);
    }
    delete gr;
}

void Tools_ScaleGraphAndMove::MoveTGraphErrorsTest() {
    const size_t n = 5;
    double x[n]  = { 1, 2, 3, 4, 5 };
    double ex[n] = { 0.1, 0.1, 0.1, 0.1, 0.1 };
	double y[n]  = { 5, 4, 3, 2, 1 };
    double ey[n] = { 0.1, 0.1, 0.1, 0.1, 0.1 };
    
    TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);
    
    MathTools::ScaleGraphAndMove(gr, 1.0, -5.0);
    
    for(uint i = 0; i<n; ++i){
        double _x, _y;
        double _ex = gr->GetErrorX(i);
        gr->GetPoint(i, _x, _y);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( x[i] - 5.0, _x, 0.00001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( ex[i], _ex, 0.00001);
    }
    delete gr;
    
}
