/*
 * sigs.h
 *
 *  Created on: Oct 12, 2014
 *      Author: hanh
 */

#ifndef SIGS_H_
#define SIGS_H_

#if defined (__cplusplus)
extern "C" {
#endif

typedef struct {
	int hashLen;
	const unsigned char *hash;

	int signatureLen;
	const unsigned char *signature;

	int pubkeyLen;
	const unsigned char *pubkey;
} Sig;

int *secp256k1_ecdsa_verify_batch(int sigsLen, const Sig *sigs);

#if defined (__cplusplus)
}
#endif

#endif /* SIGS_H_ */
