package com.example.tracker;

import android.content.Context;
import android.os.Environment;
import android.util.Log;

import com.example.tracker.model.Model;
import com.example.tracker.model.Models;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class ReadTxtFile {
    // Attribute
    private Context context;

    // Constructor
    ReadTxtFile(Context baseContext) {
        this.context = baseContext;
    }

    // READ TLE-SAT DATA FILE
    void readFile() {
        if (isExternalStorageReadable()) {
            try {
                File textfile = new File(context.getExternalFilesDir(null), "amateur.txt");
                FileInputStream fis = new FileInputStream(textfile);
                InputStreamReader isr = new InputStreamReader(fis);
                BufferedReader buff = new BufferedReader(isr);
                String line;

                ArrayList<String> arrayList = new ArrayList<>();
                Models models = Models.getInstance();
                while ((line = buff.readLine()) != null) {
                    arrayList.add(line);
                }
                fis.close();

                Model.nrOfSatellites = (arrayList.size()) / 3;
                System.out.println("Nr. of satellites: "+Model.nrOfSatellites);
                System.out.println("Test: "+arrayList.get(0));

                for (int i = 0; i < arrayList.size(); ) {
                    Model model = new Model();
                    model.setSatName(arrayList.get(i).trim());
                    model.setLineOne(arrayList.get(i + 1));
                    model.setLineTwo(arrayList.get(i + 2));
                    models.setModels(model); //In ArrayList speichern
                    i += 3;
                }
                return ;
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return ;
    }

    private boolean isExternalStorageReadable() {
        if (Environment.MEDIA_MOUNTED.equals(Environment.getExternalStorageState()) || Environment.MEDIA_MOUNTED_READ_ONLY.equals(Environment.getExternalStorageState())) {
            Log.i("State", "Is readable!");
            return true;
        } else {
            return false;
        }
    }
}
