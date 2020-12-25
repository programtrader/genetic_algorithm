#include <iostream>     //IO流
#include <cmath>        //数学类库
#include <vector>       //向量库
#include <ctime>        //时间库
#include <map>          //映射

int lower_b,upper_b;


///binary to  decimal
int decimal(std::string num){//十进制
    int decimal = 0;//十进制变量等于0
//    printf("size = %d\n",num.size());//打印一下看字符串长度
    for (int i = 0; i < num.size(); ++i) {//循环，次数为二级制数字长度。
        int power = (int) pow(2, (num.size() - i - 1) );//pow()函数用来求x的y次幂（次方）
        char c = num[i];//取出要计算的位数
        ///这里就是char c = 1,2,3,4
        int base = c - 0;//简单的转换就能变成数字
        decimal += base * power;//累加各位数
    }
//    printf("decimal = %d\n",decimal);//输出转换后的数字
    return decimal;//返回值
}

///decode chromosome//编码 染色体
double decode(int low_bound,int upper_bound,std::string chromosome){
    int deno = static_cast<int>( pow( 2 ,chromosome.size() ) ) - 1;//static_cast<int>强制转换int，2的染色体长度次方-1
    double decimal_num = static_cast<double>( decimal(chromosome) );//static_cast<double>强制转换double，decimal转换为10进制
    double x = static_cast<double>(low_bound)//low边界
              + static_cast<double>(decimal_num  * (upper_bound - low_bound) / deno );//+十进制数据*（upper边界-low边界）/整个样本空间
    return x;
}
///产生种群
void gene_chromosome(std::vector<std::string>& chromo_vec,int vec_count,int chromo_size){//产生种群
    //参数 chromo_vec 染色体向量
    //int vec_count,向量数量
    //int chromo_size 染色体的长度。
    srand((unsigned)time(NULL));//初始化，伪随机数发生器，种子为时间
    //产生一定数量的染色体
    for (int i = 0; i < vec_count; ++i) {//以向量数量为循环次数
        std::string chromosome;//染色体
        //产生固定长度的染色体
        for (int j = 0; j < chromo_size; ++j) {
            char bit = static_cast<char> (rand() % 2);//用2整除的余数编码染色体各个位数
            chromosome+=bit;//一比特信息
        }
        chromo_vec.push_back(chromosome);//将染色体推入向量中去
    }
}
///接下来的演变包括三个部分，选择(select),交叉(crossover),变异(mutation)
///对种群中的染色体进行适应度函数计算，如果是计算函数，也就是让f(x)最大的值
double f_x(std::string chromosome){//计算目标函数值
    double decoded_number = decode(lower_b,upper_b,chromosome);//先解码染色体
    double fx = decoded_number + 10*std::sin(5*decoded_number) + 7*std::cos(4*decoded_number);//然后在求值
    return fx;//返回函数值
}
///选择种群中表现好的,retain_rate，选择多大的适应性强的染色体进行保留.random_select_rate从原本应该淘汰的中选择的比例
void select(std::vector<std::string>& parents_chromo,std::vector<std::string>& chromo_vec,double remain_rate = 0.8,double random_select_rate = 0.5){
    ///传入最后一个参数让map按key降序排列
    /// 参数：
    /// parents_chromo向量双亲染色体
    /// double remain_rate 存活概率，默认0.8
    /// double random_select_rate 随机选择概率，默认0.5
    std::map< double,std::string,std::greater<double> > chromo_map;//染色体映射，排序函数std::greater<double>
    for (int i = 0; i < chromo_vec.size(); ++i) {//以染色体向量长度为个数开始循环
        chromo_map.insert(std::make_pair( f_x(chromo_vec[i]),chromo_vec[i] ) );//向map中插入染色体及其目标函数
    }
    ///只有前面的百分之remain_rate的染色体幸存下来
    int remain_length = static_cast<int>( chromo_map.size() * remain_rate );//remain_rate存活比率，被默认为0.8

    std::map< double,std::string,std::greater<double> >::iterator iter = chromo_map.begin();//遍历映射
    for (int j = 0; j < remain_length; ++j) {//循环，循环长度为remain_length
        parents_chromo.push_back(iter->second);//推入双亲队列
        ++iter;//遍历器++
    }
    for (int k = remain_length ; k < chromo_map.size(); ++k) {//循环，从remain_length,到染色体映射长度，注意++k
        if( rand() % 100 < (random_select_rate * 100) ){//产生一个随机数，让这个随指数整除100，得到的余数，如果小于自然选择概率
            parents_chromo.push_back(iter->second);//将这个随机存活的个体推入队列
        }
        ++iter;//遍历器自增
    }

}
///对上一步中幸存的parent进行交叉，产生后代
void crossover(std::vector<std::string> parents_chromo,std::vector<std::string>& chromo_vec){
    /// 函数交叉
    /// 参数：
    /// parents_chromo 亲代向量
    /// chromo_vec 种群向量
    int needed_size = chromo_vec.size() - parents_chromo.size();//需要的大小=种群向量大小-亲代向量大小
    std::vector<std::string> need_generated;//声明一个向量 需要发生
    while(need_generated.size() < needed_size){//当需要发生的大小小于需要的大小
        int father_index = rand() % parents_chromo.size();//父亲从双亲队列中随机找一个出来
        int mother_index = rand() % parents_chromo.size();//母亲从双亲队列中随机找一个出来
        //要是同一个怎么办？无性繁殖了
        std::string need_chromo;//定义一个 新的染色体个体
        if(father_index != mother_index){//如果父亲不等于母亲个体
            int random_pos = rand() % chromo_vec[0].size();//
            need_chromo.insert(
                need_chromo.end(),
                parents_chromo[father_index].begin(),
                parents_chromo[father_index].end() - random_pos);
            //插入need_chromo字符串(在指定位置loc前插入区间[start, end)的所有元素,其实只有father一个元素，在最后节点前，这样插入的是一个副本
            need_chromo.insert(need_chromo.end(),
                parents_chromo[mother_index].begin() + (parents_chromo[mother_index].size() - random_pos),
                parents_chromo[mother_index].end());
            //再插入一个mother
        }
        need_generated.push_back(need_chromo);//将need_chromo推入需要生成vector
    }
    chromo_vec.assign(parents_chromo.begin(),parents_chromo.end());//重新分配内存，先把亲代向量装进去
    chromo_vec.insert(chromo_vec.end(),need_generated.begin(),need_generated.end());//再把产生的子代向量装进去
//    printf("generated size = %d\n",chromo_vec.size());//输出子代向量的大小
    //但这个函数里拿来的交叉呢？向来和带insert的两句有关系
    //原来第二句有点长，后面包含一个随机的交叉点。
    //这是一个随机的交叉点
    //把代码重新调整一下
}
void mutation(std::vector<std::string>& chromo_vec,int mutation_rate = 0.01){
    ///函数：突变
    /// 参数：
    /// 种群向量
    /// 突变概率，默认0.01
    /// 
    for (int i = 0; i < chromo_vec.size(); ++i) {
        if( rand() % 100 < (mutation_rate * 100) ){
            int random_pos = rand() % chromo_vec[0].size();
            if( chromo_vec[i].at(random_pos) == '0')
                chromo_vec[i].at(random_pos) = '1';
            else if(chromo_vec[i].at(random_pos) == '1')
                chromo_vec[i].at(random_pos) = '0';
        }
    }
}
void evolve(std::vector<std::string>& chromo_vec){
    std::vector<std::string> parents_vec;
    select(parents_vec,chromo_vec);
    crossover(parents_vec,chromo_vec);
    mutation(chromo_vec);
}
double result(std::vector<std::string>& chromo_vec){
    std::map< double,std::string,std::greater<double> > chromo_map;
    for (int i = 0; i < chromo_vec.size(); ++i) {
        chromo_map.insert(std::make_pair( f_x(chromo_vec[i]),chromo_vec[i] ) );
    }
    printf("result = %.4f\n",decode(lower_b,upper_b,chromo_map.begin()->second));
    return decode(lower_b,upper_b,chromo_map.begin()->second);
}
int main() {
    ///例题 f(x) = x + 10*sin(5*x) + 7*cos(4*x)  在区间[0,9]的最大值，函数图像在附件
    ///[0,9] 我们假设求精度为0.001的解,则有 2^11<3000<2^12，所以每个染色体大小为12
    std::vector<std::string> chromo_vec;
    ///种群大小
    int population_count = 300;
    ///每个染色体大小
    int chromosome_size = 12;
    ///上界下界
    lower_b = 0;
    upper_b = 9;
    gene_chromosome(chromo_vec,population_count,chromosome_size);
    ///迭代次数
    int iteration_count = 100;
    for (int i = 0; i < iteration_count; ++i) {
        evolve(chromo_vec);
    }
    double res = result(chromo_vec);
    double fx = res + 10*std::sin(5*res) + 7*std::cos(4*res);
    printf("%.4f\n",fx);
    std::cout << "Hello, World!" << std::endl;
    return 0;
}