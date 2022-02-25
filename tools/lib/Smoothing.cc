#include "Smoothing.hh"

namespace SiFi
{
namespace tools
{

TH2F* SmoothGauss(TH2F* hin, double sigma)
{

    if (sigma <= 0)
    {
        std::cout << "Smearing with sigma = " << sigma
                  << " will not work, provide a positive number here..." << std::endl;
        return nullptr;
    }

    TH2F* hout = dynamic_cast<TH2F*>(hin->Clone(Form("%s_smooth", hin->GetName())));
    hout->Reset();
    const int nbinsx = hin->GetNbinsX();
    const int nbinsy = hin->GetNbinsY();
    double binwx = hin->GetXaxis()->GetBinWidth(1);
    double binwy = hin->GetYaxis()->GetBinWidth(1);

    double kernelx[nbinsx];
    double kernely[nbinsy];

    TF1* gaus = new TF1("gaus", "gaus");
    gaus->SetParameters(1. / sqrt(TMath::TwoPi()) / sigma, 0, sigma);

    for (int i = 0; i < nbinsx; i++)
        kernelx[i] = gaus->Eval(binwx * i);
    for (int i = 0; i < nbinsy; i++)
        kernely[i] = gaus->Eval(binwy * i);

    int deltabin = 0;
    double z = 0;

    // smearing in rows
    for (int biny = 1; biny < nbinsy; biny++)
    {
        for (int binx = 1; binx < nbinsx; binx++)
        {
            z = 0;
            for (int binxp = 1; binxp < nbinsx; binxp++)
            {
                deltabin = abs(binxp - binx);
                z += kernelx[deltabin] * hin->GetBinContent(binxp, biny);
            }
            hout->SetBinContent(binx, biny, z);
        }
    }
    TH2F* htmp = dynamic_cast<TH2F*>(hout->Clone());
    hout->Reset();

    // smearing in columns
    for (int binx = 1; binx < nbinsx; binx++)
    {
        for (int biny = 1; biny < nbinsy; biny++)
        {
            z = 0;
            for (int binyp = 1; binyp < nbinsy; binyp++)
            {
                deltabin = abs(binyp - biny);
                z += kernely[deltabin] * htmp->GetBinContent(binx, binyp);
            }
            hout->SetBinContent(binx, biny, z);
        }
    }

    return hout;
}

std::vector<Double_t> UQI_MSE(TH2F* sourceImage, TH2F* recoImage)
{
    std::vector<Double_t> params;
    // spdlog::info("Hello");
    // params[0] = 0.0;
    // params[1] = 0.0;
    params.push_back(0.0);
    params.push_back(0.0);
    if (sourceImage->GetNbinsX() != recoImage->GetNbinsX() ||
        sourceImage->GetNbinsY() != recoImage->GetNbinsY())
    {
        spdlog::error("Inconsistent parameters of source and reco histogram");
        // exit(EXIT_FAILURE);
        return params;
    }

    int N = recoImage->GetNbinsX() * recoImage->GetNbinsY();
    // Normalization
    double sourcemax = sourceImage->GetMaximum();
    sourceImage->Scale(1 / sourcemax);
    double recomax = recoImage->GetMaximum();
    double recomin = recoImage->GetMinimum();

    double sourceSum = 0;
    double recoSum = 0;
    double sourceVal, recoVal;

    double mse = 0;

    for (int i = 0; i < N; i++)
    {
        recoVal = (recoImage->GetBinContent(i) - recomin) / (recomax - recomin);
        sourceVal = sourceImage->GetBinContent(i);
        sourceSum += sourceVal;
        recoSum += recoVal;
        mse += pow(recoVal - sourceVal, 2);
    }
    double sourceMean = sourceSum / N;
    double recoMean = recoSum / N;

    double sourceSigma = 0;
    double recoSigma = 0;
    double cov = 0;

    for (int i = 0; i < N; i++)
    {
        sourceVal = sourceImage->GetBinContent(i);
        recoVal = (recoImage->GetBinContent(i) - recomin) / (recomax - recomin);

        sourceSigma += pow(sourceVal - sourceMean, 2);
        recoSigma += pow(recoVal - recoMean, 2);
        cov += (sourceVal - sourceMean) * (recoVal - recoMean);
    }
    sourceSigma /= (N - 1);
    recoSigma /= (N - 1);
    cov /= (N - 1);

    params[0] = mse / N;
    params[1] = 4 * cov * sourceMean * recoMean / (sourceSigma + recoSigma) /
                (pow(sourceMean, 2) + pow(recoMean, 2));
    spdlog::info("MSE = {}", params[0]);
    spdlog::info("UQI = {}", params[1]);
    return params;
}

} // namespace tools
} // namespace SiFi
