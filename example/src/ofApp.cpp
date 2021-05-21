#include "ofApp.h"

void ofApp::setup(){
	randomize();
}

void ofApp::draw(){
	ofBackground(120);

	// Draw polygons
	for (ofxPolyPol polygon : polygons) {
		polygon.draw();
	}

	// Draw lines
	ofSetColor(0);
	if (drawLines) {
		for (ofxPolyLine line : lines) {
			line.display();
		}
	}

	// Text
	ofSetColor(255);
	ofDrawBitmapString("Click to randomize", 100, 100);
	ofDrawBitmapString("Press to hide/show the lines", 100, 120);
	ofDrawBitmapString("nPolygons: " + ofToString(polygons.size()), 100, 140);
}

void ofApp::randomize() {
	lines.clear();

	for (int i = 0; i < N_LINES; i++) {
		float startX = ofRandom(ofGetWidth());
		float startY = ofRandom(ofGetHeight());
		float endX = ofRandom(ofGetWidth());
		float endY = ofRandom(ofGetHeight());

		lines.push_back({ { startX, startY }, { endX, endY } });
	}
	
	polygons = polygonDetector.detectPolygons(lines);
}

void ofApp::keyPressed(int key){
	drawLines = !drawLines;
}

void ofApp::mousePressed(int x, int y, int button){
	randomize();
}
