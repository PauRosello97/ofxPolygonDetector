#pragma once

#include "ofMain.h"
#include "ofxPolygonDetector.h";

class ofApp : public ofBaseApp{
	int N_LINES = 20;
	ofxPolygonDetector polygonDetector;
	vector<ofxPolyLine> lines;
	vector<ofxPolyPol> polygons;
	bool drawLines = true;

	public:
		void setup();
		void draw();
		void keyPressed(int key);
		void mousePressed(int x, int y, int button);
		void randomize();
};
