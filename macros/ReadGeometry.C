Bool_t ReadGeometry(TString name = "CC")
{
    gSystem->Load("libGeom");
    gSystem->Exec("pwd");
    if (name == "CC")
        TGeoManager::Import("../results/CCSimulation_TGeometry_gen4.root");
    else if (name = "CM")
        TGeoManager::Import("../results/CMsimulation_1_geometry.root");
    else
    {
        cout << name << "- incorrect option for file opening..." << endl;
        return kFALSE;
    }
    TView3D* view = (TView3D*)TView::CreateView(1);
    // gGeoManager->SetTopVisible();
    gGeoManager->GetTopVolume()->Draw();
    view->ShowAxis();
    return kTRUE;
}
