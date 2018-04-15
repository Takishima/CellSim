#include <iostream>

namespace cellsim {
namespace test {
     int test_koenigsberger(int, char**);
     int test_model(int, char**);
     int test_gap_junctions(int, char**);
     int test_small_scale(int, char**);
     int test_rk4(int, char**);
     int test_ftcs(int, char**);
     int test_model_equations(int, char**);
}
}

int main(int argc, char** argv)
{
     using namespace cellsim::test;

     // std::cout << "Testing Koenigsberger Model..." << std::endl << std::endl;

     // test_koenigsberger(argc, argv);

     // std::cout << "\n========================================\n";
     
     // std::cout << "Testing 2009 Model..." << std::endl << std::endl;

     // test_model(argc, argv);

     // std::cout << "\n========================================\n";

     // std::cout << "Testing gap junctions..." << std::endl << std::endl;

     // test_gap_junctions(argc, argv);

     // std::cout << "\n========================================\n";

     // std::cout << "Testing small scale system..." << std::endl << std::endl;

     // test_small_scale(argc, argv);

     // std::cout << "\n========================================\n";

     // std::cout << "Testing RK4..." << std::endl << std::endl;
     
     // test_rk4(argc, argv);
     
     // std::cout << "\n========================================\n";

     // std::cout << "Testing FTCS..." << std::endl << std::endl;
     
     // test_ftcs(argc, argv);
     
     // std::cout << "\n========================================\n";

     std::cout << "Testing model equations II..." << std::endl << std::endl;
     
     test_model_equations(argc, argv);
     
     std::cout << "\n========================================\n";

     return 0;
}
