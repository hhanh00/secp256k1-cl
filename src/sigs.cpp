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
#include <CL/cl.h>

#include "sigs.h"

using namespace std;

typedef unsigned char byte;


typedef struct {
	string hash;
	string signature;
	string pubkey;
} Sigcpp;

class Verifier {
private:
	cl_context context;
	cl_command_queue command_queue;
	cl_program program;
public:
	Verifier(cl_context context, cl_command_queue command_queue, cl_program program) :
		context(context), command_queue(command_queue), program(program) {}
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

	cl_int ret;
    cl_kernel precomp_kernel = clCreateKernel(program, "secp256k1_ecmult_table_precomp_gej", &ret);
    cl_kernel mult_kernel = clCreateKernel(program, "secp256k1_ecmult", &ret);
	// Everything must check ...
	int *rets = secp256k1_ecdsa_verify_batch(context, command_queue, precomp_kernel, mult_kernel, sigs.size(), &*sigs.begin());
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

    FILE *fp;
    char *source_str;
    size_t source_size;

    // Get platform and device information
    cl_platform_id platform_id = NULL;
    cl_device_id device_id = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1,
            &device_id, &ret_num_devices);

    // Create an OpenCL context
    cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

    cl_program program;
    {
		ifstream source("secp256k1.cl", ios_base::in|ios_base::binary|ios::ate);
		size_t size = source.tellg();
		auto source_str = boost::make_unique<char []>(size);
		const char *sources[] = { source_str.get() };
		const size_t lens[] = { size };
		source.seekg(0, ios::beg);
		source.read(source_str.get(), size);

	    program = clCreateProgramWithSource(context, 1,
	    		sources, lens, &ret);

	    // Build the program
	    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
	    cout << ret << endl;

	    size_t len;
		char *buffer;
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
		buffer = (char *)malloc(len);
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
		cout << buffer << endl;
    }

	ifstream sigs("/tmp/sigs.dat", ios_base::in|ios_base::binary|ios::ate);
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
			if (i == 262144)
				break;
		}
	}
	cout << toVerify.size() << endl;

	Verifier v(context, command_queue, program);

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

	// CPU: cout << elapsed << endl; // 190s for 2142967 sigs => 0.08ms / sig


	return 0;
}
