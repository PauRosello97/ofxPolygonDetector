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
		void update();
		void draw();
		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		void randomize();
};
