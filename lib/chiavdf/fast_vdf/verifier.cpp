#include <flint/fmpz.h>
#include "include.h"
#include "integer_common.h"
#include "vdf_new.h"
#include "picosha2.h"
#include "nucomp.h"
#include "proof_common.h"
#include "create_discriminant.h"

form FastPowFormNucomp(form &x, integer &D, integer num_iterations, integer &L)
{
    if (num_iterations == integer(0))
        return form::identity(D);

    integer new_num_iterations = num_iterations;
    new_num_iterations >>= 1;
    form res = FastPowFormNucomp(x, D, new_num_iterations, L);
    nucomp_form(res, res, res, D, L);
    if (num_iterations % integer(2) == integer(1))
        nucomp_form(res, res, x, D, L);
    res.reduce();
    return res;
}

integer ConvertBytesToInt(uint8_t *bytes, int start_index, int end_index)
{
    integer res(0);
    bool negative = false;
    if (bytes[start_index] & (1 << 7))
        negative = true;
    for (int i = start_index; i < end_index; i++)
    {
        res = res * integer(256);
        if (!negative)
            res = res + integer(bytes[i]);
        else
            res = res + integer(bytes[i] ^ 255);
    }
    if (negative)
    {
        res = res + integer(1);
        res = res * integer(-1);
    }
    return res;
}

form DeserializeForm(integer &d, uint8_t *bytes, int int_size)
{
    integer a = ConvertBytesToInt(bytes, 0, int_size);
    integer b = ConvertBytesToInt(bytes, int_size, 2 * int_size);    
    form f = form::from_abd(a, b, d);
    return f;
}

std::vector<form> DeserializeProof(uint8_t *proof_bytes, int proof_len, integer &D)
{
    int int_size = (D.num_bits() + 16) >> 4;
    std::vector<form> proof;
    for (int i = 0; i < proof_len; i += 2 * int_size)
    {
        std::vector<uint8_t> tmp_bytes;
        for (int j = 0; j < 2 * int_size; j++)
            tmp_bytes.push_back(proof_bytes[i + j]);
        proof.emplace_back(DeserializeForm(D, tmp_bytes.data(), int_size));
    }
    return proof;
}

void VerifyWesolowskiProof(integer &D, form x, form y, form proof, int iters, bool &is_valid)
{
    int int_size = (D.num_bits() + 16) >> 4;
    integer L = root(-D, 4);
    integer B = GetB(D, x, y);
    integer r = FastPow(2, iters, B);
    form f1 = FastPowFormNucomp(proof, D, B, L);
    form f2 = FastPowFormNucomp(x, D, r, L);
    if (f1 * f2 == y)
    {
        is_valid = true;
    }
    else
    {
        is_valid = false;
    }
}

bool CheckProofOfTimeNWesolowskiInner(integer &D, form x, uint8_t *proof_blob,
                                      int blob_len, int iters, int int_size,
                                      std::vector<int> iter_list, int recursion)
{
    uint8_t result_bytes[10000];
    uint8_t proof_bytes[10000];
    memcpy(result_bytes, proof_blob, 2 * int_size);
    memcpy(proof_bytes, proof_blob + 2 * int_size, blob_len - 2 * int_size);
    form y = DeserializeForm(D, result_bytes, int_size);
    std::vector<form> proof = DeserializeProof(proof_bytes, blob_len - 2 * int_size, D);
    if (recursion * 2 + 1 != proof.size())
        return false;
    if (proof.size() == 1)
    {
        bool is_valid;
        VerifyWesolowskiProof(D, x, y, proof[0], iters, is_valid);
        return is_valid;
    }
    else
    {
        if (!(proof.size() % 2 == 1 && proof.size() > 2))
            return false;
        int iters1 = iter_list[iter_list.size() - 1];
        int iters2 = iters - iters1;
        bool ver_outer;
        std::thread t(VerifyWesolowskiProof, std::ref(D), x, proof[proof.size() - 2], proof[proof.size() - 1], iters1, std::ref(ver_outer));
        uint8_t new_proof_bytes[10000];
        for (int i = 0; i < blob_len - 4 * int_size; i++)
            new_proof_bytes[i] = proof_blob[i];
        iter_list.pop_back();
        bool ver_inner = CheckProofOfTimeNWesolowskiInner(D, proof[proof.size() - 2], new_proof_bytes, blob_len - 4 * int_size, iters2, int_size, iter_list, recursion - 1);
        t.join();
        if (ver_inner && ver_outer)
            return true;
        return false;
    }
}

bool CheckProofOfTimeNWesolowski(integer &D, form x, uint8_t *proof_blob, int proof_blob_len, int iters, int recursion)
{
    int int_size = (D.num_bits() + 16) >> 4;
    uint8_t new_proof_blob[10000];
    int new_cnt = 4 * int_size;
    memcpy(new_proof_blob, proof_blob, new_cnt);
    std::vector<int> iter_list;
    for (int i = new_cnt; i < proof_blob_len; i += 4 * int_size + 8)
    {
        auto iter_vector = ConvertBytesToInt(proof_blob, i, i + 8).to_vector();
        iter_list.push_back(iter_vector[0]);
        if (iter_vector[0] < 0)
            return false;
        memcpy(new_proof_blob + new_cnt, proof_blob + i + 8, 4 * int_size);
        new_cnt += 4 * int_size;
    }
    bool is_valid = CheckProofOfTimeNWesolowskiInner(D, x, new_proof_blob, new_cnt, iters, int_size, iter_list, recursion);
    return is_valid;
}

std::vector<uint8_t> HexToBytes(char *hex_proof)
{
    int len = strlen(hex_proof);
    assert(len % 2 == 0);
    std::vector<uint8_t> result;
    for (int i = 0; i < len; i += 2)
    {
        int hex1 = hex_proof[i] >= 'a' ? (hex_proof[i] - 'a' + 10) : (hex_proof[i] - '0');
        int hex2 = hex_proof[i + 1] >= 'a' ? (hex_proof[i + 1] - 'a' + 10) : (hex_proof[i + 1] - '0');
        result.push_back(hex1 * 16 + hex2);
    }
    return result;
}

// Intended to match with ProofOfTime type from Chia Blockchain.
struct ProofOfTimeType
{
    integer discriminant;
    integer a;
    integer b;
    uint64_t iterations_needed;
    std::vector<uint8_t> witness;
    uint8_t witness_type;

    ProofOfTimeType(const integer &discriminant, const integer &a, const integer &b, uint64_t iterations_needed,
                    const std::vector<uint8_t> &witness, uint8_t witness_type)
    {
        this->discriminant = discriminant;
        this->a = a;
        this->b = b;
        this->iterations_needed = iterations_needed;
        this->witness = witness;
        this->witness_type = witness_type;
    }
};

// Converts from ProofOfTimeType to CheckProofOfTimeNWesolowski-like arguments and calls the check.
bool CheckProofOfTimeType(ProofOfTimeType &proof)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    bool result;

    try
    {
        form x = form::generator(proof.discriminant);
        int int_size = (proof.discriminant.num_bits() + 16) >> 4;
        form y = form::from_abd(proof.a, proof.b, proof.discriminant);
        std::vector<uint8_t> proof_blob = SerializeForm(y, int_size);
        proof_blob.insert(proof_blob.end(), proof.witness.begin(), proof.witness.end());

        result = CheckProofOfTimeNWesolowski(proof.discriminant, x, proof_blob.data(), proof_blob.size(), proof.iterations_needed, proof.witness_type);
    }
    catch (std::exception &e)
    {
        result = false;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    print("Verification time = " + to_string(duration) + "ms.");

    return result;
}

int main()
{
    auto t1_start = std::chrono::high_resolution_clock::now();
    char challenge_hash1_hex[] = "a4bb1461ade74ac602e9ae511af68bb254dfe65d61b7faf9fab82d0b4364a30b";
    auto challenge_hash1 = HexToBytes(challenge_hash1_hex);
    integer discriminant1 = CreateDiscriminant(challenge_hash1);
    auto t1_end = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t1_end - t1_start).count();
    print("Create discriminant time = " + to_string(duration1) + "ms.");
    print("Discriminant = ");
    print(discriminant1.impl);

    auto t2_start = std::chrono::high_resolution_clock::now();
    char challenge_hash2_hex[] = "1633f29c0ca0597258507bc7d323a8bd485d5f059da56340a2c616081fb05b7f";
    auto challenge_hash2 = HexToBytes(challenge_hash2_hex);
    integer discriminant2 = CreateDiscriminant(challenge_hash2);
    auto t2_end = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2_end - t2_start).count();
    print("Create discriminant time = " + to_string(duration2) + "ms.");
    print("Discriminant = ");
    print(discriminant2.impl);

    // Test 1: block 11278 challenge_hash = a4bb1461ade74ac602e9ae511af68bb254dfe65d61b7faf9fab82d0b4364a30b
    std::vector<uint8_t> witness1(
        {0, 66, 83, 222, 194, 35, 255, 62, 48, 148, 1, 100, 235, 199, 156, 156,
         24, 228, 16, 253, 25, 189, 42, 83, 211, 151, 30, 214, 205, 234, 176,
         143, 103, 11, 236, 108, 48, 90, 73, 253, 252, 123, 98, 49, 100, 211,
         134, 126, 122, 253, 197, 214, 156, 104, 29, 235, 77, 255, 37, 110, 150,
         132, 26, 76, 14, 255, 190, 212, 244, 102, 121, 29, 91, 97, 152, 214,
         228, 231, 118, 31, 184, 33, 250, 213, 245, 240, 182, 126, 20, 92, 240,
         178, 227, 83, 39, 7, 84, 122, 138, 82, 188, 99, 80, 63, 141, 215, 155,
         45, 140, 155, 191, 236, 204, 62, 110, 69, 240, 196, 139, 117, 240, 69,
         153, 205, 15, 11, 78, 73, 66, 113, 0, 0, 0, 0, 0, 19, 108, 68, 0, 18,
         56, 13, 141, 27, 60, 129, 110, 67, 132, 125, 163, 147, 28, 132, 246, 13,
         5, 246, 1, 214, 193, 6, 6, 83, 5, 31, 41, 40, 23, 100, 134, 110, 234,
         191, 175, 164, 89, 152, 158, 210, 170, 191, 181, 231, 141, 115, 187,
         188, 117, 13, 16, 118, 57, 14, 200, 9, 38, 155, 56, 207, 29, 199, 184,
         255, 242, 183, 152, 131, 47, 139, 66, 200, 145, 3, 230, 113, 118, 112,
         96, 31, 147, 55, 96, 0, 139, 108, 111, 245, 44, 88, 171, 220, 228, 40,
         255, 245, 85, 110, 163, 132, 156, 190, 217, 61, 224, 83, 141, 162, 183,
         185, 176, 34, 220, 62, 140, 161, 193, 130, 254, 192, 33, 171, 38, 190,
         151, 23, 62, 45, 0, 26, 202, 156, 87, 237, 11, 119, 153, 67, 48, 84,
         126, 74, 210, 174, 20, 26, 114, 242, 211, 158, 56, 88, 206, 124, 247,
         157, 227, 206, 218, 111, 157, 21, 109, 45, 216, 194, 12, 192, 72, 81,
         141, 99, 120, 135, 212, 244, 37, 97, 53, 59, 205, 178, 153, 244, 36,
         173, 68, 79, 191, 22, 150, 104, 243, 0, 6, 67, 238, 203, 215, 17, 96,
         97, 51, 19, 126, 9, 150, 201, 139, 128, 132, 38, 119, 124, 215, 139,
         217, 4, 125, 75, 52, 216, 180, 80, 26, 47, 153, 151, 98, 153, 35, 16,
         158, 184, 185, 84, 136, 248, 61, 227, 105, 149, 46, 157, 207, 36, 132,
         128, 10, 103, 246, 199, 197, 156, 10, 197, 193, 163, 0, 0, 0, 0, 0, 58,
         69, 48, 0, 101, 133, 189, 101, 125, 9, 91, 174, 190, 12, 26, 159, 234,
         224, 17, 36, 112, 170, 14, 206, 164, 160, 20, 140, 144, 250, 67, 81,
         231, 68, 65, 172, 145, 188, 239, 49, 78, 48, 178, 167, 87, 102, 14, 21,
         183, 126, 141, 86, 143, 75, 163, 175, 202, 96, 7, 177, 176, 112, 239,
         41, 178, 222, 118, 225, 255, 219, 166, 43, 222, 69, 28, 31, 16, 187,
         213, 112, 113, 190, 240, 227, 141, 175, 195, 218, 52, 95, 236, 58, 122,
         173, 23, 142, 171, 222, 33, 155, 232, 18, 247, 119, 139, 51, 218, 202,
         37, 181, 50, 60, 78, 214, 164, 89, 244, 30, 190, 11, 115, 56, 153, 170,
         154, 239, 139, 143, 50, 100, 239, 85, 141, 0, 33, 14, 241, 123, 70, 18,
         20, 190, 103, 31, 183, 124, 146, 5, 254, 69, 120, 191, 173, 56, 219,
         126, 111, 177, 223, 13, 181, 75, 155, 121, 84, 163, 110, 185, 145, 245,
         8, 191, 32, 241, 228, 98, 162, 77, 4, 139, 199, 24, 135, 4, 75, 165, 9,
         2, 101, 117, 49, 83, 13, 233, 160, 90, 54, 48, 0, 18, 224, 163, 145, 90,
         8, 46, 176, 190, 175, 151, 24, 190, 8, 238, 94, 208, 74, 164, 186, 90,
         17, 164, 243, 22, 151, 77, 36, 61, 7, 15, 215, 4, 65, 187, 12, 134, 51,
         91, 114, 229, 146, 13, 219, 233, 99, 27, 168, 167, 117, 179, 115, 147,
         119, 118, 154, 147, 143, 186, 148, 218, 207, 108, 221});
    // Test 2: block 11279 challenge_hash = 1633f29c0ca0597258507bc7d323a8bd485d5f059da56340a2c616081fb05b7f
    std::vector<uint8_t> witness2(
        {0, 95, 166, 136, 156, 117, 76, 69, 226, 152, 64, 70, 228, 158, 116, 8,
         12, 143, 234, 213, 47, 134, 1, 130, 242, 28, 56, 196, 193, 38, 249, 40,
         116, 127, 167, 205, 140, 210, 50, 216, 126, 195, 87, 55, 16, 139, 169,
         98, 92, 160, 54, 55, 99, 49, 203, 2, 52, 204, 136, 45, 54, 163, 147, 99,
         192, 255, 168, 3, 227, 201, 122, 245, 198, 189, 183, 188, 218, 16, 61,
         177, 46, 216, 197, 59, 245, 120, 250, 29, 173, 203, 65, 116, 18, 41, 92,
         67, 35, 91, 143, 220, 36, 62, 90, 185, 57, 159, 126, 27, 232, 33, 75,
         227, 141, 170, 67, 37, 1, 243, 127, 66, 204, 246, 252, 252, 65, 97, 178,
         7, 221, 147, 0, 0, 0, 0, 0, 52, 24, 184, 0, 52, 80, 239, 39, 35, 177,
         86, 60, 172, 107, 216, 35, 74, 230, 230, 33, 50, 81, 14, 170, 245, 99,
         198, 31, 152, 230, 12, 230, 217, 74, 165, 227, 154, 137, 190, 211, 212,
         28, 138, 149, 225, 108, 15, 130, 199, 153, 150, 49, 121, 26, 247, 239,
         28, 160, 139, 135, 115, 201, 33, 24, 116, 103, 80, 252, 255, 214, 129,
         147, 33, 40, 197, 86, 41, 10, 157, 158, 114, 86, 155, 235, 245, 232,
         131, 147, 153, 0, 103, 94, 168, 217, 236, 175, 190, 37, 171, 120, 200,
         162, 31, 28, 114, 162, 58, 242, 82, 105, 68, 125, 7, 163, 86, 144, 5,
         95, 197, 243, 72, 134, 24, 66, 2, 253, 161, 164, 254, 8, 112, 10, 91, 0,
         41, 133, 59, 14, 93, 97, 223, 50, 96, 61, 100, 48, 13, 165, 183, 134,
         42, 170, 193, 99, 33, 181, 76, 17, 130, 26, 100, 145, 110, 220, 236, 79,
         253, 182, 213, 165, 168, 131, 112, 49, 85, 116, 226, 248, 95, 8, 172,
         79, 154, 26, 89, 100, 216, 105, 89, 57, 70, 215, 147, 248, 176, 164,
         102, 253, 0, 4, 170, 116, 152, 110, 74, 25, 185, 120, 146, 206, 191,
         178, 51, 180, 134, 19, 199, 131, 100, 88, 54, 183, 45, 40, 236, 58, 197,
         51, 23, 92, 83, 204, 202, 105, 78, 181, 153, 8, 220, 68, 15, 47, 250,
         136, 170, 224, 232, 86, 170, 169, 254, 145, 172, 139, 62, 240, 174, 6,
         68, 28, 100, 219, 175, 0, 0, 0, 0, 0, 156, 74, 40, 0, 68, 114, 66, 105,
         155, 227, 40, 194, 10, 114, 70, 170, 51, 16, 127, 59, 111, 184, 157,
         228, 172, 245, 145, 226, 113, 220, 60, 186, 58, 57, 113, 149, 133, 21,
         157, 71, 116, 40, 139, 4, 138, 70, 10, 27, 93, 151, 118, 126, 104, 22,
         136, 187, 156, 191, 237, 77, 203, 114, 184, 85, 139, 132, 115, 16, 0, 3,
         179, 45, 126, 100, 104, 90, 115, 91, 68, 225, 195, 129, 169, 168, 58,
         89, 251, 244, 106, 162, 182, 154, 209, 54, 194, 101, 46, 255, 65, 137,
         243, 80, 10, 155, 205, 186, 195, 143, 125, 116, 255, 186, 191, 102, 242,
         50, 5, 225, 137, 97, 219, 232, 58, 184, 246, 54, 164, 50, 252, 126, 50,
         10, 173, 0, 70, 78, 200, 97, 172, 104, 190, 76, 17, 170, 127, 90, 13,
         207, 65, 153, 125, 126, 119, 246, 155, 176, 212, 219, 79, 16, 175, 38,
         176, 82, 88, 92, 254, 143, 236, 11, 37, 61, 155, 99, 126, 200, 10, 51,
         64, 178, 179, 152, 35, 250, 171, 117, 79, 180, 214, 115, 193, 50, 57,
         197, 125, 174, 65, 254, 255, 224, 183, 19, 125, 0, 80, 163, 94, 29, 30,
         206, 94, 136, 112, 151, 78, 91, 236, 159, 171, 57, 112, 137, 2, 0, 195,
         140, 48, 186, 193, 211, 121, 191, 22, 222, 172, 126, 50, 194, 80, 93,
         69, 13, 47, 79, 46, 159, 226, 19, 154, 49, 212, 159, 83, 230, 248, 200,
         32, 37, 189, 169, 219, 245, 43});

    ProofOfTimeType test1(
        ///*discriminant=*/integer("-146034004988995324450899632927708569076216311660629486868969524905823791101129965616306403106228613155248048777324859112377567377611365244344698891812959891333662090324001598409511351166803295083408004513522410571274209999557757565527745109011222985028131502058812870261063827579868894671524386074886059562807"),
        /*discriminant=*/discriminant1,
        /*a=*/integer("2120439809548002603699476515242856052401016261690016272875375071691961393091525387070479560771020084644071207726495978942882631782725802522199105524842492"),
        /*b=*/integer("-1447745607279856836390430406078195096118390105960950865534703611396791372115445411292284152604500691817421284831093402993306669441174757567119801738100995"),
        /*iterations_needed=*/5728275,
        /*witness=*/witness1,
        /*witness_type=*/2
    );

    ProofOfTimeType test2(
        ///*discriminant=*/integer("-141140317794792668862943332656856519378482291428727287413318722089216448567155737094768903643716404517549715385664163360316296284155310058980984373770517398492951860161717960368874227473669336541818575166839209228684755811071416376384551902149780184532086881683576071479646499601330824259260645952517205526679"),
        /*discriminant=*/discriminant2,
        /*a=*/integer("5995822174752598681464821178400044775263491831121191947378966754053931067055682917999005598919776725834807951347990196033583183503215350056672265407022667"),
        /*b=*/integer("1857599280246213629471341186238029632043329881356694152473161696212923061132534499554981549765761965469743565522045552935302619566637106283783998935421377"),
        /*iterations_needed=*/15363962,
        /*witness=*/witness2,
        /*witness_type=*/2
    );

    bool result1 = CheckProofOfTimeType(test1);
    if (result1)
        print("Test #1 passed.\n");
    else
        print("Test #1 failed.\n");
    
    bool result2 = CheckProofOfTimeType(test2);
    if (result2)
        print("Test #2 passed.\n");
    else
        print("Test #2 failed.\n");

    // Test 3, modified correct proof by 1 byte.

    witness1[witness1.size() - 1]++;
    ProofOfTimeType test3(
        /*discriminant=*/integer("-146034004988995324450899632927708569076216311660629486868969524905823791101129965616306403106228613155248048777324859112377567377611365244344698891812959891333662090324001598409511351166803295083408004513522410571274209999557757565527745109011222985028131502058812870261063827579868894671524386074886059562807"),
        /*a=*/integer("2120439809548002603699476515242856052401016261690016272875375071691961393091525387070479560771020084644071207726495978942882631782725802522199105524842492"),
        /*b=*/integer("-1447745607279856836390430406078195096118390105960950865534703611396791372115445411292284152604500691817421284831093402993306669441174757567119801738100995"),
        /*iterations_needed=*/5728275,
        /*witness=*/witness1,
        /*witness_type=*/2
    );
    bool result3 = CheckProofOfTimeType(test3);
    if (!result3)
        print("Test #3 passed.\n");
    else
        print("Test #3 failed.\n");
    
    // Test 4, modified iterations needed.

    ProofOfTimeType test4(
        /*discriminant=*/integer("-141140317794792668862943332656856519378482291428727287413318722089216448567155737094768903643716404517549715385664163360316296284155310058980984373770517398492951860161717960368874227473669336541818575166839209228684755811071416376384551902149780184532086881683576071479646499601330824259260645952517205526679"),
        /*a=*/integer("5995822174752598681464821178400044775263491831121191947378966754053931067055682917999005598919776725834807951347990196033583183503215350056672265407022667"),
        /*b=*/integer("1857599280246213629471341186238029632043329881356694152473161696212923061132534499554981549765761965469743565522045552935302619566637106283783998935421377"),
        /*iterations_needed=*/15363963,
        /*witness=*/witness2,
        /*witness_type=*/2
    );

    bool result4 = CheckProofOfTimeType(test4);
    if (!result4)
        print("Test #4 passed.\n");
    else
        print("Test #4 failed.\n");

    return 0;
}