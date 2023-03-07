// HomeWork2 FFT Iterative and Parallel
// Name: Chaitanya Attarde SUID: 660307939

#include <iostream>
#include <complex>
#include <math.h>
#include <vector>
#include <fstream>
#include <thread>
#include <chrono>

using namespace std;
using namespace chrono;
using namespace literals::chrono_literals;
typedef complex<double> cx;

class Timer {
public:

    chrono::system_clock::time_point Begin, End;
    chrono::system_clock::duration Runtime;
    Timer() {
        Begin = chrono::system_clock::now();
    }
    ~Timer() {
        End = chrono::system_clock::now();
        Runtime = End - Begin;
        cout << "Runtime = " << chrono::duration_cast<chrono::milliseconds>(Runtime).count() << " ms" << endl;
    }
};

const double PI = 3.1415926536;
unsigned int bitReverse(unsigned int x, int n);
void fft_Itr(vector<cx>& a,int n1,int n2);
void fft_common(vector<cx>& a , int n1 , int n2, int inc_size , int part_size);
void fft_recursive_parallel(vector<cx> &a,int n1,int n2);

int main() {
    vector<cx> a;     //for recursive FFT
    vector<cx> a1;    //for iterative FFT
    
    // Read Data
    ifstream In("input1024.txt");
    if (!In) {
        cout << "Fail to open file" << endl;
        return 0;
    }
    double d1;
    while (In >> d1) {
        a.emplace_back(cx(d1));
        a1.emplace_back(cx(d1));
    }
    In.close();
    
    //bit reversal
    for(int i=0; i < int(a.size()); i++){
        int x = bitReverse(i, int(a.size()));
        if(x > i){
            swap(a[i],a[x]);
            swap(a1[i],a1[x]);
        }
    }
    
    //Recursive Parallel FFT
    {
        Timer TT;
        int n1 = 0,n2= int(a.size())-1;
        thread t1(fft_recursive_parallel,ref(a),n1,(n1+n2)/2);
        fft_recursive_parallel(a,(n1+n2+1)/2,n2);
        if(t1.joinable())t1.join();
        
        // threaded for last_step
        {
            thread t1(fft_common, ref(a) , 0 , int(a.size())-1 , 1, 0);
            thread t2(fft_common, ref(a), 0 , int(a.size())-1,  1 , 1);
            thread t3(fft_common, ref(a), 0 , int(a.size())-1, 1 , 2);
            fft_common(a , 0 , int(a.size())-1 , 1 , 3);
            if(t1.joinable())t1.join();
            if(t2.joinable())t2.join();
            if(t3.joinable())t3.join();
        }
    }
    
    // Iterative Parallel FFT
    {
        Timer TT1;
        int n = int(a.size());
        // Run for input less than or equal to 8
        if(n <= 8){
            fft_Itr(a1,0,n-1);
        }
        
        // Run paralley for input greater than 8
        else{
            // 4 threads for FFT : 0  to log2(n)-2
            {
                thread t1(fft_Itr,ref(a1),0,(n/4)-1);
                thread t2(fft_Itr,ref(a1),n/4,(n/2)-1);
                thread t3(fft_Itr,ref(a1),n/2, (3*(n/4))-1);
                fft_Itr(a1,3*(n/4),n-1);
                if(t1.joinable())t1.join();
                if(t2.joinable())t2.join();
                if(t3.joinable())t3.join();
            }
            // 4 threads for FFT : 0  to log2(n)-1
            {
                thread t1(fft_common, ref(a1) , 0 , (n/2) - 1 , 2 , 0 );
                thread t2(fft_common, ref(a1), 0 , (n/2) - 1 , 2 , 1 );
                thread t3(fft_common, ref(a1), n/2 , n-1 , 2 ,  0 );
                fft_common(a1 , n/2 , n-1 , 2 , 1 );
                if(t1.joinable())t1.join();
                if(t2.joinable())t2.join();
                if(t3.joinable())t3.join();
            }
            // 4 threads for FFT : Final Step
            {
                thread t1(fft_common, ref(a1) , 0 , n-1 , 1, 0);
                thread t2(fft_common, ref(a1), 0 , n-1,  1 , 1);
                thread t3(fft_common, ref(a1), 0 , n-1, 1 , 2);
                fft_common(a1 , 0 , n-1 , 1 , 3);
                if(t1.joinable())t1.join();
                if(t2.joinable())t2.join();
                if(t3.joinable())t3.join();
            }
        }
    }
    
    ofstream Out("output_data_recursive.txt");//Use the right file name described in HW2
    ofstream Out1("output_data_iterative.txt");//Use the right file name described in HW2
    if (!Out || !Out1) {
        cout << "Fail to access file ...." << endl;
        return 0;
    }
    for (int i = 0; i < a.size(); ++i){
        Out << a[i] << "\n";
        Out1 << a1[i] << "\n";
    }
    Out.close();
    return 0;
}

unsigned int bitReverse(unsigned int num, int n)
{
    int logn = 0;
    string str1 = "";
    while(logn <= (log2(n) - 1)){
        str1 += ((num >> logn) & 0x1) == 0 ? '0' : '1';
        logn++;
    }
    return stoi(str1,nullptr,2);

}

void fft_recursive_parallel(vector<cx> &a,int n1,int n2){
    if((n2 - n1) == 0)return;
    if((n2 - n1 + 1) == int(a.size()/2)){
        thread t2(fft_recursive_parallel,ref(a),n1,(n1+n2)/2);
        fft_recursive_parallel(a,(n1+n2+1)/2,n2);
        if(t2.joinable())t2.join();
        // threaded for last second_last_step
        thread t1(fft_common, ref(a) , n1 , n2 , 2, 0);
        fft_common(a,n1,n2,2,1);
        if(t1.joinable())t1.join();
    }
    else{
        fft_recursive_parallel(a,n1,(n1+n2)/2);
        fft_recursive_parallel(a,(n1+n2+1)/2,n2);
        for(int i=0; i<=(n2-n1)/2;i++){
            int com_size = (n2-n1 +1)/2;
            int inc_size = int(a.size()/(n2-n1+1));
            cx temp = a[n1+i],J{0,1};
            cx WW = i==0?1:exp(cx(-2)*J*cx(PI)*cx(inc_size * i)/cx(int(a.size())));
            a[i+ n1] = temp + (WW * a[n1+ i + com_size]);
            a[n1+ i + com_size] = temp - (WW * a[n1+ i + com_size]);
        }
    }
    
}

void fft_Itr(vector<cx>&a, int n1,int n2){
    int N = n2-n1+1;
    for(int i=2; i<=(n2-n1+1); i=i*2){
        int j = 0, inc_size = int(a.size()/i);
        while(j < N){
            for(int k = 0; k < (i/2); k++){
                cx temp = a[ n1 + k + j ],J{0,1};
                cx WW = k==0?1:exp(cx(-2)*J*cx(PI)*cx(inc_size * k)/cx(int(a.size())));
                a[ n1 + k + j ] = temp + (WW * a[ n1 + k + j + i/2]);
                a[ n1 + k + j + i/2] = temp - (WW * a[ n1 + k + j + i/2]);
            }
            j+=i;
        }
    }
}

void fft_common(vector<cx> &a, int n1 , int n2 , int inc_size , int part_size){
    int j = (part_size * int(a.size()/8));
    int com_size = (n2 - n1 + 1)/2;
    for(int i = j ; i < (j + int(a.size()/8)) ;i++){
        cx temp = a[i + n1],J{0,1};
        cx WW = i==0?1:exp(cx(-2)*J*cx(PI)*cx(inc_size * i)/cx(a.size()));
        a[i + n1] = temp + (WW * a[i + n1 + com_size]);
        a[i + n1 + com_size] = temp - (WW * a[i + n1 + com_size]);
    }
}



