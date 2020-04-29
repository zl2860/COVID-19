library(tidyverse)
library(numDeriv)
library(pracma)

# Import & Clean data
####################################################################################
data = read_csv("./covid19-1.csv") %>% janitor::clean_names()
#skimr::skim(data)

#看一眼国家的数据
#names(data)
#length(unique(data$country_region)) #169
#data %>%
#  group_by(country_region) %>%
#  summarise(n_cases = max(confirmed_cases)) %>% #confirmed_cases = 累计cases
#  filter(n_cases != 0) %>%
#  arrange(desc(n_cases)) %>% #nrow(.) = 145
#  .[1:10,] %>%
#  ggplot(aes(x = reorder(country_region,-n_cases), y = log(n_cases),fill = #country_region)) +
#  geom_col() +
#  theme_bw() +
#  ggsci::scale_fill_d3() +
#  labs(title = "Top 10 Countries With the Most Cases") +
#  theme(plot.title = element_text(hjust = .5))#  145个国家有cases,top 10

#有cases的国家汇总
data_havecases = data %>%
  group_by(country_region) %>%
  summarise(n_cases = max(confirmed_cases)) %>% #confirmed_cases = 累计cases
  filter(n_cases != 0)

country = unique(data_havecases$country_region) # 挑出来有cases 的国家

#整理日期数据，准备拿来优化
data = data %>% 
  mutate(date = str_c(date,"20"), 
         date = as.Date.character(date,
                                  format = "%m/%d/%Y")) %>%
  filter(country_region %in% country) #这个data 改过日期，剔除了没有cases的国家


#############################################################
# 定义一个函数：先选定一个国家，提取数据，然后把日期转换为天数
date_to_day = function(df,country_name){
  df_country = df %>%
    filter(country_region == country_name) %>% #选出这个地区的数据
    group_by(date) %>% 
    summarise(cases = sum(confirmed_cases)) %>% # 按日期求累积cases的和，得出每天整个国家的累积cases
    arrange(date) #按日期排序
  
  start_day = which(df_country$cases!=0)[1] # 第一次出现感染的日期对应的行数
  start_date = df_country$date[start_day] # 第一次感染的日期
  df_country = df_country[start_day:nrow(df_country),] %>%
    mutate(days_since_spread = c(1:nrow(.))) # 创建logistic curve 里的t
  
  return(list(start_date = start_date,
              df_country = df_country,
              country.name = country_name))  
} # 返回spread date & 从那一天起感染的数据 & 国家名（用df_country来优化）
#############################################################


# example
t = date_to_day(data,"Portugal")
t$start_date
t$df_country
t$country.name

# 打印10个国家的数据
#for (region in country[1:10]) {
#  res = date_to_day(data,region)
#  print(res$df_country)
#}

t

plot(t$df_country$days_since_spread,t$df_country$cases)
####################################################################################



# 定义优化用到的函数
####################################################################################
# log.curve
#创建一个logistic curve
log.curve = function(para,t){
  a = para[1]
  b = para[2]
  c = para[3]
  return(a/(1+exp(-b*(t - c))))
}

# N_i+1 - N_i 的向量,计算每个时间点cases差值
delta = function(x){
  r = NULL
  for (i in 2:length(x)) {
    r[i-1] = x[i] - x[i-1]
  }
  return(r)
}

# 检查delta，原始数据包含连续4天或者以上没有增长的caes的时候返回True
check_delta = function(delta){
  if(sum(min(delta[1:3])) == 0){
    return(TRUE) #as.integer((2/5)*length(delta))
  }else{
    return(FALSE)
  }
}

#先变个型 a <= k, bc <= a, b <= r
#先估计k,数据太垃圾的时候（连续好多天没增长）就选择用拐点法，否则用三点法
k.estimate = function(df){ # para = c(k,a,r)
  cases = df$cases
  n1 = cases[1]
  n2 = cases[as.integer(median(1:length(cases)))] #2 * t2 = t3?
  n3 = cases[length(cases)]
  #print(c(n1,n2,n3))
  delta = delta(cases)
  #print(delta)
  if(check_delta(delta) == FALSE & ((2*n1*n2*n3)-n2^2*(n1+n3))/(n1*n3-n2^2) > 0){
    return(((2*n1*n2*n3)-n2^2*(n1+n3))/(n1*n3-n2^2))
  }
    #vec = NULL
    ##if(n3 == max(cases)){return(2*n3)}
    ##else{
    #for (i in 2:length(cases)) {
    #  vec[i-1] = cases[i] - cases[i-1]
    #}
    if(cases[which.max(delta)+1] != max(cases) && abs(cases[which.max(delta)+1] - max(cases)) > (1/7) * cases[which.max(delta)+1]){
      return(2*max(cases))
    }else{
    return(2*cases[which.max(delta)+1])
      }#+2??
  }

#根据估计的k所得出真实Y
y.real = function(df){
  N = df$cases
  k = rep(k.estimate(df), length(N))
  for (i in 1:length(k)) {
    if (k[i] == N[i]) {
      k[i] = k[i] + 0.5
    }
  }
  trans = (k - N)/N
  y.real = NULL
  
  for (each in trans) {
    if (each > 0){
      each = log(each)
      y.real = append(y.real,each)
    }else{
      each = 1/2 * log(each^2)
      y.real = append(y.real,each)
    }
  }
  return(y.real)
}

#拟合值Y-hat
y.hat = function(df,para){ # c(a,r)
  a = para[1]
  r = para[2]
  t = df$days_since_spread
  return(a - r*t)
}


# loss function
loss = function(y.real,y.hat,para){
  return((1/length(y.real)) * sum((y.real - y.hat)^2))
}

# 求梯度
gradient = function(para,y.real,y.hat){ #para = c(a,r)
  grad_a = (-2/length(y.real)) * sum(y.real - y.hat)
  grad_r = (-2/length(y.real)) * sum((y.real - y.hat)*-c(1:length(y.real))) # 此处可能会有错
  return(c(grad_a,grad_r)) # 返回grad
}


#检查数据
check_df = function(df){
  if(sum(df$cases == min(df$cases))>=10|sum(df$cases == min(df$cases)+1)>=10|sum(df$cases == min(df$cases)+2) >= 10){
    return(TRUE) #垃圾数据，需要裁剪
  }else{
    return(FALSE)
  }
}

#定义一个函数，先检测数据的分布，然后决定是否裁剪。如果早期有很多没有增长的点，则进行裁剪再丢进optimizer，否则可以直接丢进去
df.modifier = function(df){
  if(check_df(df) == TRUE){
    if(sum(df$cases == min(df$cases))>=10){s = sum(df$cases == min(df$cases))}
    if(sum(df$cases == min(df$cases)+1)>=10){s = sum(df$cases == min(df$cases)+1)}
    if(sum(df$cases == min(df$cases)+2)>=10){s = sum(df$cases == min(df$cases)+2)}
    slicing_index = s - 3
    df = df %>% .[-c(1:slicing_index),] %>% mutate(days_since_spread = 1:nrow(.))
    #重新改日子了免得影响拟合
    return(df)
  }else{
    return(df)
  }
}





# 优化过程


#优化器，梯度下降优化a,r; 返回值为老师给的公式里的a,b,c,loss,iter
optimizer = function(para, df, tol = 1e-10, maxiter = 20000){
  k = k.estimate(df)
  y.real = y.real(df)
  y.hat = y.hat(df,para)
  
  loss_prev = -Inf
  loss_cur = loss(y.real,y.hat,para) # loss
  para_cur = para
  iter = 1
  res = data.frame(a = para_cur[1], r = para_cur[2], loss = loss_cur,iter = iter)
  
  while (abs(loss_cur - loss_prev) > tol && iter < maxiter) {
    iter = iter + 1
    #print(iter)
    step = 1
    #print(para_cur)
    loss_prev = loss_cur
    para_prev = para_cur
    
    grad = gradient(para_prev, y.real = y.real, y.hat = y.hat)
    
    #print(grad)
    
    para_cur = para_prev - step * grad # update parameters
    
    #print(para_cur)
    
    y.hat = y.hat(df,para_cur) #更新参数后的y.hat
    #print(y.hat)
    
    loss_cur = loss(y.real = y.real, y.hat = y.hat,para_cur)#更新参数后的loss
    
    #print(paste("loss_cur > loss_pre:", loss_cur > loss_prev))
    #print(loss_cur)
    while (loss_cur > loss_prev) {
      step = 0.5 * step
      para_cur = para_prev - step * grad
      y.hat = y.hat(df,para_cur)
      loss_cur = loss(y.real = y.real, y.hat = y.hat,para_cur)
    }
    
    res_tmp  = data.frame(a = para_cur[1], r = para_cur[2], loss = loss_cur,iter = iter)
    res = rbind(res,res_tmp)
    
  }
  res = res %>% mutate(k = k,b = r,c = a/r) %>% select(loss,k,b,c,iter) %>% rename("a" = "k")
  return(res)
}
####################################################################################



# obtain parameters
####################################################################################
##获取单个国家abc
country #所有可用国家的列表
obtain_para = function(country){
  t = date_to_day(data,country)
  df = df.modifier(t$df_country)
  para = c(1,4) # 随缘设定的a,r起始值，别改了就这样吧
  
  para.data = optimizer(para,df) 
  para.final = para.data[nrow(para.data),c(2:4)] %>% as.numeric() #最后的参数
  pred = log.curve(para = para.final, t = 1:nrow(df))
  return(list(para = para.final, pred = pred, time = 1:nrow(df), turncated.df = df))
}
####################################################################################

# (a,b,c) for all countries
####################################################################################
para.all = function(country){
  para.all = data.frame()
  for (country.each in country) {
    
    a = obtain_para(country.each)$para[1]
    b = obtain_para(country.each)$para[2]
    c = obtain_para(country.each)$para[3]
    country.name = country.each
    para.all = rbind(para.all, data.frame(a = a, b = b, c = c, country = country.name))
    print(paste(country.each," done!"))
  }
  return(para.all)
}

#test = country[1:5] # 取前5个国家
#parameters = para.all(test) #结果
####################################################################################


# 画图测试
####################################################################################
#  随便pick一个国家画图，目前只看到泰国不是很ok，原因是泰国的数据属于轻微垃圾水平，之前定义的检测垃圾数据的函数没有鉴别出来所有没有进行裁剪。数据太不规则的国家先不管。
country.plot = "US"#sample(country,1)
res.test = obtain_para(country.plot)
t = date_to_day(data,country.plot)
df = t$df_country
print(country.plot)
plot(res.test$turncated.df$cases)
lines(res.test$time,res.test$pred)
#k.estimate(df)
####################################################################################


df = date_to_day(data,"Antigua and Barbuda")$df_country
df$cases
df = df.modifier(df)
k.estimate(df)
country = country[-113]
parameters = para.all(country)
write.csv(parameters,"./parameters_0427.csv")
