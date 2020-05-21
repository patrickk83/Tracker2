package com.example.tracker.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class Models {

    private static Models instance;
    private ArrayList<Model> models = new ArrayList<>();

    private Models () {}

    public static Models getInstance () {
        if (Models.instance == null) {
            Models.instance = new Models ();
        }
        return Models.instance;
    }

       public ArrayList<Model> getModels() {
        return models;
    }

    public void setModels(Model model) {
        models.add(model);
    }

    public  void sortToName(){
        Collections.sort(models, Comparator.comparing(Model::getSatName));
    }

}


