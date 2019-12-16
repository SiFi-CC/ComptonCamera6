#include <cppunit/extensions/HelperMacros.h>

#include <MathTools.hh>

class Tools_SmoothGauss : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( Tools_SmoothGauss );
    CPPUNIT_TEST( SmoothGaussTest );
    CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp();
    void tearDown();
    
protected:
    void SmoothGaussTest();
};

//------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( Tools_SmoothGauss );

//------------------------------------------------------------------
void Tools_SmoothGauss::setUp() {
}

//------------------------------------------------------------------
void Tools_SmoothGauss::tearDown() {
}

//------------------------------------------------------------------

void Tools_SmoothGauss::SmoothGaussTest() {
    
    TH2F *hOrig = new TH2F("hOrig", "hOrig", 100, -50, 50, 100, -50, 50);
    int binCenter = hOrig->FindBin(0,0);
    hOrig->SetBinContent(binCenter, 1000);
    
    double sigma = 10;
    
    TH2F *hSmear = dynamic_cast<TH2F*>(MathTools::SmoothGauss(hOrig, sigma));
    
    TH1F *projX = (TH1F*)hSmear->ProjectionX();
    TH1F *projY = (TH1F*)hSmear->ProjectionY();
    
    TF1 *funX = new TF1("funX", "gaus", -50, 50);
    TF1 *funY = new TF1("funY", "gaus", -50, 50);

    projX->Fit(funX, "Q0");
    projY->Fit(funY, "Q0");

    CPPUNIT_ASSERT_DOUBLES_EQUAL( sigma, funX->GetParameter(2), 0.1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( sigma, funY->GetParameter(2), 0.1);
    
    delete hOrig, hSmear;
    delete projX, projY;
    delete funX, funY;
}
