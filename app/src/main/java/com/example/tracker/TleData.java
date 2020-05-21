package com.example.tracker;

public class TleData {

    //Basic attributes
    private String satName;
    private String satLineOneLineNumber;
    private String satLineOneNumber;
    private String satLineOneClass;
    private String satLineOneLaunchyr;
    private String satLineOneLaunchpc;
    private String satLineOneEpochyr;
    private String satLineOneEpochday;
    private String satLineOneFtdmm;
    private String satLineOneStdmm;
    private String satLineOneDrag;
    private String satLineOneEph;
    private String satLineOneEle;
    private String satLineOneChksum1;

    private String satLineTwoLineNumber;
    private String satLineTwoNumber;
    private String satLineTwoIncl;
    private String satLineTwoRa;
    private String satLineTwoEcc;
    private String satLineTwoPeri;
    private String satLineTwoMa;
    private String satLineTwoMm;
    private String satLineTwoRevnr;
    private String satLineTwoChecksum2;


    //Constructor
    TleData(String satName, String satLineOneLineNumber, String satLineOneNumber, String satLineOneClass, String satLineOneLaunchyr, String satLineOneLaunchpc,
            String satLineOneEpochyr, String satLineOneEpochday, String satLineOneFtdmm, String satLineOneStdmm, String satLineOneDrag, String satLineOneEph,
            String satLineOneEle, String satLineOneChksum1, String satLineTwoLineNumber, String satLineTwoNumber, String satLineTwoIncl, String satLineTwoRa,
            String satLineTwoEcc, String satLineTwoPeri, String satLineTwoMa, String satLineTwoMm, String satLineTwoRevnr, String satLineTwoChecksum2){

        this.satName = satName;
        this.satLineOneLineNumber = satLineOneLineNumber;
        this.satLineOneNumber = satLineOneNumber;
        this.satLineOneClass = satLineOneClass;
        this.satLineOneLaunchyr = satLineOneLaunchyr;
        this.satLineOneLaunchpc = satLineOneLaunchpc;
        this.satLineOneEpochyr = satLineOneEpochyr;
        this.satLineOneEpochday = satLineOneEpochday;
        this.satLineOneFtdmm = satLineOneFtdmm;
        this.satLineOneStdmm = satLineOneStdmm;
        this.satLineOneDrag = satLineOneDrag;
        this.satLineOneEph = satLineOneEph;
        this.satLineOneEle = satLineOneEle;
        this.satLineOneChksum1 = satLineOneChksum1;
        this.satLineTwoLineNumber = satLineTwoLineNumber;
        this.satLineTwoNumber = satLineTwoNumber;
        this.satLineTwoIncl = satLineTwoIncl;
        this.satLineTwoRa = satLineTwoRa;
        this.satLineTwoEcc = satLineTwoEcc;
        this.satLineTwoPeri = satLineTwoPeri;
        this.satLineTwoMa = satLineTwoMa;
        this.satLineTwoMm = satLineTwoMm;
        this.satLineTwoRevnr = satLineTwoRevnr;
        this.satLineTwoChecksum2 = satLineTwoChecksum2;
    }

    //Getters
    public String getSatName() {return satName;}
    public String getSatLineOneLineNumber() {return satLineOneLineNumber;}
    public String getSatLineOneNumber() {return satLineOneNumber;}
    public String getSatLineOneClass() {return satLineOneClass;}
    public String getSatLineOneLaunchyr() {return satLineOneLaunchyr;}
    public String getSatLineOneLaunchpc() {return satLineOneLaunchpc;}
    public String getSatLineOneEpochyr() {return satLineOneEpochyr;}
    public String getSatLineOneEpochday() {return satLineOneEpochday;}
    public String getSatLineOneFtdmm() {return satLineOneFtdmm;}
    public String getSatLineOneStdmm() {return satLineOneStdmm;}
    public String getSatLineOneDrag() {return satLineOneDrag;}
    public String getSatLineOneEph() {return satLineOneEph;}
    public String getSatLineOneEle() {return satLineOneEle;}
    public String getSatLineOneChksum1() {return satLineOneChksum1;}
    public String getSatLineTwoLineNumber() {return satLineTwoLineNumber;}
    public String getSatLineTwoNumber() {return satLineTwoNumber;}
    public String getSatLineTwoIncl() {return satLineTwoIncl;}
    public String getSatLineTwoRa() {return satLineTwoRa;}
    public String getSatLineTwoEcc() {return satLineTwoEcc;}
    public String getSatLineTwoPeri() {return satLineTwoPeri;}
    public String getSatLineTwoMa() {return satLineTwoMa;}
    public String getSatLineTwoMm() {return satLineTwoMm;}
    public String getSatLineTwoRevnr() {return satLineTwoRevnr;}
    public String getSatLineTwoChecksum2() {return satLineTwoChecksum2;}

}
