#include <cppunit/extensions/HelperMacros.h>

#include <MathTools.hh>

class Tools_SubtractBackground : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( Tools_SubtractBackground );
    CPPUNIT_TEST( SubtrBackground_TF1 );
    CPPUNIT_TEST( SubtrBackground_Spline );
    CPPUNIT_TEST( SubtrBackground_Graph );
	CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp();
    void tearDown();
    
protected:
    void SubtrBackground_TF1();
    void SubtrBackground_Spline();
    void SubtrBackground_Graph();
    
private:
    TGraphErrors *fGraph;
    TGraphErrors *fGraphBg;
};

//------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( Tools_SubtractBackground );

//------------------------------------------------------------------
void Tools_SubtractBackground::setUp() {
    
    const int npoints = 20;
    
    //----- creating test graph
    double x[npoints]  =  { -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, 
                            -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5 };
    double y[npoints]  =  { 3, 3, 3, 3, 3, 3, 4, 23, 144, 334, 357, 128, 27, 6, 4, 3, 3, 3, 3, 3 };
    double ex[npoints] =  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    double ey[npoints] =  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    
    fGraph = new TGraphErrors(npoints, x, y, ex, ey);

    //----- creating background graph
    double _x[npoints]  =  { -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, 
                             -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5 };
    double _y[npoints]  =  { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };
    double _ex[npoints] =  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    double _ey[npoints] =  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    
    fGraphBg = new TGraphErrors(npoints, _x, _y, _ex, _ey);
}

//------------------------------------------------------------------
void Tools_SubtractBackground::tearDown() {
}

//------------------------------------------------------------------
void Tools_SubtractBackground::SubtrBackground_TF1() {
    
    const int npoints = 20;
    
    TF1 *bgfun = new TF1("bgfun", "pol1", -10, 10);
    bgfun->FixParameter(0, 3.0834);
    bgfun->FixParameter(1, 7.2334E-3);
    
    TGraphErrors *graphSig = MathTools::SubtractBackground(fGraph, bgfun);
    
    for(int i=0; i<npoints; i++){
        
        double sigX, sigY;
        double origX, origY;
        
        graphSig->GetPoint(i, sigX, sigY);
        fGraph->GetPoint(i, origX, origY);
        
        CPPUNIT_ASSERT_DOUBLES_EQUAL( bgfun->Eval(sigX)+sigY, origY, 0.00001 );
    }
    
    delete graphSig, bgfun;
}

//------------------------------------------------------------------
void Tools_SubtractBackground::SubtrBackground_Spline() {
    
    const int npoints = 20;
    
    TSpline3 *spline = new TSpline3("spline", fGraphBg);
    TGraphErrors *graphSig = MathTools::SubtractBackground(fGraph, spline);
    
    for(int i=0; i<npoints; i++){
        
        double sigX, sigY;
        double origX, origY;
        
        graphSig->GetPoint(i, sigX, sigY);
        fGraph->GetPoint(i, origX, origY);
        
        CPPUNIT_ASSERT_DOUBLES_EQUAL( spline->Eval(sigX)+sigY, origY, 0.00001 );
    }
    
    delete graphSig, spline;
}

//------------------------------------------------------------------
void Tools_SubtractBackground::SubtrBackground_Graph() {
    
    const int npoints = 20;
    
    TGraphErrors *graphSig = MathTools::SubtractBackground(fGraph, fGraphBg);
    
    for(int i=0; i<npoints; i++){
        
        double sigX, sigY;
        double origX, origY;
        
        graphSig->GetPoint(i, sigX, sigY);
        fGraph->GetPoint(i, origX, origY);
        
        CPPUNIT_ASSERT_DOUBLES_EQUAL( fGraphBg->Eval(sigX)+sigY, origY, 0.00001 );
    }
    
    delete graphSig;
}
