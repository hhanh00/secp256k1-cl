//============================================================================
// Name        : sigs.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cassert>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <boost/make_unique.hpp>
#include <secp256k1.h>
#include <windows.h>

#include "sigs.h"

using namespace std;

typedef unsigned char byte;


typedef struct {
	string hash;
	string signature;
	string pubkey;
} Sigcpp;

class Verifier {
public:
	bool verify(const Sig &sig);
	void verifyBatch(const vector<Sigcpp> &sigs);
};

bool Verifier::verify(const Sig &sig) {
	return secp256k1_ecdsa_verify(
			(byte *)sig.hash, sig.hashLen,
			(byte *)sig.signature, sig.signatureLen,
			(byte *)sig.pubkey, sig.pubkeyLen) == 1;
}

void Verifier::verifyBatch(const vector<Sigcpp> &sigcpps) {
	vector<Sig> sigs;
	for (auto &s : sigcpps) {
		sigs.push_back(Sig { s.hash.size(), (byte*)s.hash.data(), s.signature.size(), (byte*)s.signature.data(), s.pubkey.size(), (byte*)s.pubkey.data() });
	}
	// Everything must check
	int *rets = secp256k1_ecdsa_verify_batch(sigs.size(), &*sigs.begin());
	for (int i = 0; i < sigcpps.size(); i++) {
		if (!rets[i])
			cout << "Failed #" << i << endl; // #5 should fail
	}

	free(rets);
}


string readString(char *&p) {
	char sz = *p++;
	string s(p, p + sz);
	p += sz;
	return s;
}

int main() {
	secp256k1_start();
	ifstream sigs("/tmp/sigs2.dat", ios_base::in|ios_base::binary|ios::ate);
	auto size = sigs.tellg();
	vector<Sigcpp> toVerify;
	{
		auto memblock = boost::make_unique<char []>(size);
		sigs.seekg(0, ios::beg);
		sigs.read(memblock.get(), size);
		sigs.close();

		int i = 0;
		auto p = memblock.get();
		while (p != memblock.get() + size) {
			string hash = readString(p);
			string signature = readString(p);
			string key = readString(p);
			if (i == 5)
				hash.at(0) = 4;
			Sigcpp sig { hash, signature, key };

			toVerify.push_back(sig);
			i++;
			if (i == 10)
				break;
		}
	}
	cout << toVerify.size() << endl;

	Verifier v;

	LARGE_INTEGER q1, q2, ticksPerSecond;
	QueryPerformanceFrequency(&ticksPerSecond);
	QueryPerformanceCounter(&q1);

	v.verifyBatch(toVerify);

#if 0
	int i = 0;
	for (auto &s : toVerify) {
		auto check = v.verify(s);
		if (!check)
			cout << "Failed #" << i << endl;
		i++;
	}
#endif

	QueryPerformanceCounter(&q2);
	auto elapsed = (q2.QuadPart - q1.QuadPart) / ticksPerSecond.QuadPart;
	cout << elapsed << endl; // 190s for 2142967 sigs => 0.08ms / sig
	return 0;
}
