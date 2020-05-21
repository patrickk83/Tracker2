package com.example.tracker;

import androidx.annotation.NonNull;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.Manifest;
import android.annotation.SuppressLint;
import android.content.Context;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.location.Location;
import android.location.LocationManager;
import android.os.Bundle;
import android.provider.Settings;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;
import android.widget.Toast;

import com.example.tracker.helper.DownloadHelper;
import com.example.tracker.interfaces.DownloadInterface;
import com.example.tracker.model.Model;
import com.example.tracker.model.Models;
import com.google.android.gms.location.FusedLocationProviderClient;
import com.google.android.gms.location.LocationCallback;
import com.google.android.gms.location.LocationRequest;
import com.google.android.gms.location.LocationResult;
import com.google.android.gms.location.LocationServices;
import com.google.android.gms.tasks.OnCompleteListener;
import com.google.android.gms.tasks.Task;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;


/***************************************************************************************************
 * TEST - GPS + Download Text-File + Read Text File + Save in ArrayList
 **************************************************************************************************/

public class MainActivity extends AppCompatActivity {
    private static final String TAG = "Tracker";
    static boolean checkFirstTime = true; // Check ob Datei bereits existiert

    private static final int PERMISSION_ID = 44; //notwendig für GPS
    FusedLocationProviderClient mFusedLocationClient;

    String strLocator = "";
    CoordToLoc coordToLoc;

    double lat, lon, alt;
    float fLat, fLon;

    TextView tvMaidenhead;
    TextView tvLat;
    TextView tvLon;
    TextView tvAlt;
    TextView tvAlist;
    Models models;


    private ArrayList<TleData> tleDataArrayList = new ArrayList<>();

    @SuppressLint("SetTextI18n")
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        Button btOne = findViewById(R.id.bt_One);
        tvMaidenhead = findViewById(R.id.tv_Maid);
        tvLat = findViewById(R.id.tv_Lat);
        tvLon = findViewById(R.id.tv_Lon);
        tvAlt = findViewById(R.id.tv_Alt);
        tvAlist = findViewById(R.id.tv_ArrayList);


        //Get GPS-Coordinates
        mFusedLocationClient = LocationServices.getFusedLocationProviderClient(this);

        coordToLoc = new CoordToLoc();

        // Check if file exists
        if (checkFirstTime) {
            Log.i(TAG, "Datei wird neu geladen");
            final File file = new File(getExternalFilesDir(null), "amateur.txt");
            DownloadHelper downloadHelper = new DownloadHelper(new DownloadInterface() {
                @Override
                public void onResult() {
                    Log.i(TAG, "Datei fertig geladen und Daten werden gelesen");
                    models = Models.getInstance();
                    new ReadTxtFile(getBaseContext()).readFile();
                    readDaten();
                    models.sortToName();

                }
            }, file);
            downloadHelper.execute();
            checkFirstTime = false;

        } else {
            Log.i(TAG, "Datei ist noch actuell ");
            models = Models.getInstance();
            new ReadTxtFile(getBaseContext()).readFile();
            readDaten();
            models.sortToName();
        }
        // Button jump to Activity 2
        btOne.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                Intent intent = new Intent(MainActivity.this, HamSatActivity.class);
                startActivity(intent);
            }
        });
    }

    private void readDaten() {
        Log.i(TAG, "Daten werden gelesen");
        for (Model m : models.getModels()) {
            tleDataArrayList.add(new TleData(m.getSatName(), m.getSatLineOneLineNumber(), m.getSatLineOneNumber(), m.getSatLineOneClass(), m.getSatLineOneLaunchyr(), m.getSatLineOneLaunchpc(),
                    m.getSatLineOneEpochyr(), m.getSatLineOneEpochday(), m.getSatLineOneFtdmm(), m.getSatLineOneStdmm(), m.getSatLineOneDrag(), m.getSatLineOneEph(), m.getSatLineOneEle(),
                    m.getSatLineOneChksum1(), m.getSatLineTwoLineNumber(), m.getSatLineTwoNumber(), m.getSatLineTwoIncl(), m.getSatLineTwoRa(), m.getSatLineTwoEcc(), m.getSatLineTwoPeri(),
                    m.getSatLineTwoMa(), m.getSatLineTwoMm(), m.getSatLineTwoRevnr(), m.getSatLineTwoChksum2()));

            /*tleDataArrayList.add(new TleData("satName", "satLineOneLineNr", "satLineOneNr", "satLineOneClass",
                    "satLineOneLaunchyr", "satLineOneLaunchpc", "satLineOneEpochyr", "satLineOneEpochday",
                    "satLineOneFtdmm", "satLineOneStdmm", "satLineOneDrag", "satLineOneEph", "satLineOneEle",
                    "satLineOneChksum1", "satLineTwoLineNumber", "satLineTwoNumber", "satLineTwoIncl",
                    "satLineTwoRa", "satLineTwoEcc", "satLineTwoPeri", "satLineTwoMa", "satLineTwoMm",
                    "satLineTwoRevnr", "satLineTwoChecksum2"));
            */
        }

        //Collections.sort(tleDataArrayList, new SortById());
        //System.out.println("ArrayList[0]: " + m.get);
        //System.out.println("ArrayList[1]: "+tleDataArrayList.get(1));
        //System.out.println("ArrayList[2]: "+tleDataArrayList.get(2));
        //tvAlist.setText("ArrayList: "+tleDataArrayList);

        Collections.sort(tleDataArrayList, Comparator.comparing(TleData::getSatName));
    }


    /***********************************************************************************************
     * Methoden
     **********************************************************************************************/
    public void getLastLocation() {
        if (checkPermissions()) {
            if (isLocationEnabled()) {
                mFusedLocationClient.getLastLocation().addOnCompleteListener(
                        new OnCompleteListener<Location>() {
                            @Override
                            public void onComplete(@NonNull Task<Location> task) {
                                Location location = task.getResult();
                                if (location == null) {
                                    requestNewLocationData();
                                } else {

                                    lat = location.getLatitude();
                                    lon = location.getLongitude();
                                    alt = location.getAltitude();

                                    tvLat.setText("Latitude: " + lat);
                                    tvLon.setText("Longitude: " + lon);
                                    tvAlt.setText("Altitude: " + alt);

                                    fLat = (float) lat;
                                    fLon = (float) lon;

                                    strLocator = CoordToLoc.latLonToGridSquare(fLat, fLon);
                                    tvMaidenhead.setText(String.format("QTH: %s", strLocator));
                                }
                            }
                        }
                );
            } else {
                Toast.makeText(this, "Turn on location", Toast.LENGTH_LONG).show();
                Intent intent = new Intent(Settings.ACTION_LOCATION_SOURCE_SETTINGS);
                startActivity(intent);
            }
        } else {
            requestPermissions();
        }
    }

    private void requestNewLocationData() {
        LocationRequest mLocationRequest = new LocationRequest();
        mLocationRequest.setPriority(LocationRequest.PRIORITY_HIGH_ACCURACY);
        mLocationRequest.setInterval(0);
        mLocationRequest.setFastestInterval(0);
        mLocationRequest.setNumUpdates(1);
        mFusedLocationClient = LocationServices.getFusedLocationProviderClient(this);
        mFusedLocationClient.requestLocationUpdates(mLocationRequest, mLocationCallback, android.os.Looper.myLooper());

    }

    private LocationCallback mLocationCallback = new LocationCallback() {
        @Override
        public void onLocationResult(LocationResult locationResult) {
            Location mLastLocation = locationResult.getLastLocation();
            lat = mLastLocation.getLatitude();
            lon = mLastLocation.getLongitude();
        }
    };

    private boolean checkPermissions() {
        return ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_COARSE_LOCATION) == PackageManager.PERMISSION_GRANTED ||
                ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_FINE_LOCATION) == PackageManager.PERMISSION_GRANTED;
    }

    private void requestPermissions() {
        ActivityCompat.requestPermissions(this, new String[]{Manifest.permission.ACCESS_FINE_LOCATION, Manifest.permission.ACCESS_COARSE_LOCATION}, PERMISSION_ID);
    }

    private boolean isLocationEnabled() {
        LocationManager locationManager = (LocationManager) getSystemService(Context.LOCATION_SERVICE);
        return locationManager.isProviderEnabled(LocationManager.GPS_PROVIDER) || locationManager.isProviderEnabled(
                LocationManager.NETWORK_PROVIDER
        );
    }

    @Override
    public void onResume() {
        super.onResume();
        if (checkPermissions()) {
            getLastLocation();
        } else {
            requestPermissions();
            getLastLocation();
        }
    }

    @Override
    public void onRequestPermissionsResult(int requestCode, @NonNull String[] permissions, @NonNull int[] grantResults) {
        super.onRequestPermissionsResult(requestCode, permissions, grantResults);
        if (requestCode == PERMISSION_ID) {
            if (grantResults.length > 0 && grantResults[0] == PackageManager.PERMISSION_GRANTED) {
                getLastLocation();
            }
        }
    }


}
