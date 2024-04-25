##########################################################
## Copyright (c) NWAFU Wheat Bioincloud.lab 2022-2025
##      Project: GWAS结果批量整理
##  Description: 本脚本适用于从GWAS结果中进行显著SNP位点的整理
##         Date: 20240412
##       Author: Jewin ( zaojewin@icloud.com )
##      Version: 2.0.0
##  更新：由于基因型的格式不同，专门编写的统计程序
##########################################################
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(stringr)

# 参数设置-----

prefix <- "20240412_Test_GWAS"
windows_near <- 300*1000 #默认300kb内为连续显著区段
opts_over_maker <- 10 # 超过该值的点
# filename <- phe.model.pvalue.suffix # data文件夹内命名方式
# 注意：要求所有的文件命名方式如上

id_list <- list.files("./data/",pattern = "*.txt")
Ref <- fread("./Ref.csv") # 参考基因组信息


# 生成参数列表
id_df <- str_split(id_list,"[.]") %>% as.data.frame() %>% t() %>% as.data.frame()
type_list <- id_df$V1[!duplicated(id_df$V1)]

phe_list <- id_df$V1[!duplicated(id_df$V1)]
model_list <- c("MLM","FarmCPU")

p_value <- 6 #显著性阈值
suffix <- ".txt"

write.csv(id_df,str_c(prefix,"_parment.csv"),row.names = F)

# 创建输出文件夹
out_dirs <- c("1_GWAS_Result_txt2csv",
              "2_SNP_Infomation",
              "3_Gene_Maping_Result",
              "4_Rebind_All_Output")
for (mydir in out_dirs){
    if (dir.exists(paste0("./out/",mydir))){
        cat(paste0("[check] ./out/",mydir," >>> \t is exist !\n"))
    }else{
        dir.create(paste0("./out/",mydir))
        cat(paste0("[check] ./out/",mydir," \t create finished !\n"))
    }
}

# 转换染色体编号-用于统计染色体格式
i <- 1
new <- data.frame(matrix(ncol = 2))
new <- new[-1,]
for (x in c(1:7)){
    for (y in c("A","B","D")){
        chr <- paste0(x,y)
        new_add <- c(i,chr)
        new <- rbind(new,new_add)
        i <- i + 1
    }
}
colnames(new) <- c("CHROM","chr")

# step1：txt2csv----
for (id in id_list){
    file_name <- paste0("./data/",id)
    atom <- str_split(id,"[.]")
    type <- atom[[1]][1] %>% str_split("_") %>% .[[1]] %>% .[1]
    phe <- atom[[1]][1]
    way <- atom[[1]][2] %>% str_replace("Farm","") # 将模型替换为CPU
    plast <- p_value

    cli::cli_alert_success(str_c("Running Now is:\t",file_name,"\t",Sys.time()))
    # 计算p值并筛选
    df <- read_delim(file_name,delim = " ",skip = 1,col_names = c("INDEX","SNP","CHROM","POS","REF","ALT","Effect","SE","pvalue"),
                     col_types = cols(CHROM = col_character()))
    colnames(df)[9] <- way
    df %>% 
        mutate(log = round(-log10(!!sym(way)),1)) %>% 
        filter(log > plast) ->data
    
    # 替换染色体编号
    data2 <- data
    data2$loc <- phe
    # 待标注的log值筛选
    data2$logwt <- ifelse(data2$log > opts_over_maker,paste0('log=',data2$log,sep=""),NA)
    data2$MB <- data2$POS/1000000
    # 写出为中间结果，M表示苗期，DS、IT表示成株期
    write_csv(data2,paste0("./out/1_GWAS_Result_txt2csv/",phe,".",
                           way,".csv"))
}

# step2:comp_snp4ini----
Ref_chr <- Ref
all_single <- list() #汇总单标记
all_near <- list() #汇总连续标记
id_list_step2 <- list.files("./out/1_GWAS_Result_txt2csv/")
for (id in id_list_step2){
    # 创建染色体
    chr_list <- list()
    for (tmp_chr in new$chr){
        chr_list[[tmp_chr]] <- filter(Ref_chr,chr == tmp_chr)
    }
    cli::cli_alert_success(str_c("[Runing Model]：",id,"\t",Sys.time()))
    # 开始计算----
    file_name <- paste0("./out/1_GWAS_Result_txt2csv/",id)
    atom <- str_split(id,"[.]")
    # print(file_name)
    data <- read_csv(file_name,show_col_types = FALSE)
    loc <- data$loc[1] # 表型小种或者生态区,phe
    job <- paste0(atom[[1]][1],"_",atom[[1]][2])
    way <- atom[[1]][2]
    
    ### 单标记筛选 ========================================================================
    # 计算基因位置间距
    data$longH <- NA
    data$longQ <- NA
    data$class <- NA
    # 显著位点小于3个的情况下跳过
    if (nrow(data) < 3){
        next
    }
    for (i in 2:(nrow(data)-1)){
        a <- data$POS[i]
        i <- i+1
        b <- data$POS[i]
        i <- i-2
        c <- data$POS[i]
        i <- i+1
        
        # 更新判断SNP首尾位置的方法
        if (data$CHROM[i-1] != data$CHROM[i] &
            data$CHROM[i+1] == data$CHROM[i]){
            data$class[i] <- "shou"
            next
        }else{
            if (data$CHROM[i-1] == data$CHROM[i] &
                data$CHROM[i+1] != data$CHROM[i]){
                data$class[i] <- "wei"
                next
            }else{
                if (data$CHROM[i-1] != data$CHROM[i] &
                    data$CHROM[i+1] != data$CHROM[i]){
                    data$class[i] <- "single"
                    next
                }
            }
        }
        
        # if(i == nrow(data)){
        #     data$class[i] <- "wei"
        #     break
        # }
        # if(a-c<0){
        #     data$class[i] <- "shou"
        #     next
        # }
        # if(b-a<0){
        #     data$class[i] <- "wei"
        #     next
        # }
        data$longH[i] <- (b-a)
        data$longQ[i] <- (a-c)
    }
    data$class[1] <- "shou"
    data$class[nrow(data)] <- "wei"
    # 对距离进行区分,按照windows_near为区分阈值
    
    for (i in 1:nrow(data)){
        if (is.na(data$longH[i]) | is.na(data$longQ[i])){
            next
        }
        if (data$longH[i]>windows_near & data$longQ[i]>windows_near){
            data$class[i] <- "single"
        }
    }
    
    
    # 单标记位点处理
    data$ws <- ifelse(is.na(data$logwt),
                      paste0(data$SNP,",Find in ",str_replace(id,".csv",""),sep=""),
                      paste0(data$SNP,",Find in ",str_replace(id,".csv",""),",",data$logwt,sep=""))
    ### 单标记信息位置注释 ===================================================================
    # 单标记位置信息写入single
    single <- data.frame(matrix(ncol = 4))
    single <- single[-1,]
    colnames(single) <- c("positon","info","chr","loc") 
    for (i in 1:nrow(data)){
        tem_class <- data$class[i]
        tem_add <- c(data$POS[i],data$ws[i],data$CHROM[i],data$loc[i])
        ifelse(tem_class == 'single',single <- rbind(single,tem_add),"1")
        ifelse(tem_class == 'shou',single <- rbind(single,tem_add),"2")
        ifelse(tem_class == 'wei',single <- rbind(single,tem_add),"3")
    }
    colnames(single) <- c("positon","info","chr","loc") 
    
    ### 连续区间筛选 ===================================================================
    
    near <- data.frame(matrix(ncol = 5)) #初始化矩阵
    near <- near[-1,]
    colnames(near) <- c("p1","p2","info","chr","number") 
    
    for (x in c(1:7)){
        for (y in c("A","B","D")){
            chr_id <- paste0(x,y,sep="")
            foot <- 0 # 步长，用于迭代计算阅读框
            for (i in which(data$CHROM==chr_id)){
                if (sum(data$CHROM==chr_id)<2){next} #若某个染色体的位点数小于3则跳过
                if (foot>0){ # 如果foot变量大于0，说明两个位点存在跨越关系，进行归零
                    foot <- foot - 1 # 
                    next
                }
                n_pos_1 <- data$POS[i]
                for (m in which(data$CHROM==chr_id)){
                    n_pos_2 <- data$POS[m]
                    if (n_pos_2 - n_pos_1 < 0){ # 后一个值小于前一个值时跳过
                        next
                    }else{
                        if (n_pos_2 - n_pos_1 == 0){ # 两个值相等时跳过
                            next
                        }else{
                            if (n_pos_2 - n_pos_1 < windows_near){ # 任意两个位点距离小于预设窗口大小
                                foot <- foot+1 # 向前进行一步，跨越一个位点
                                if (is.na(data$class[i])){
                                    data$class[i] <- "near" # 如果此时位点尚不属于single、shou、wei，则为near
                                }
                            }else{
                                if (foot > 10){
                                    n_wn <- paste0(data$SNP[i],"-",data$SNP[m-1],", [",foot,"],Find in ",job)
                                }else{
                                    n_wn <- paste0(data$SNP[i],"-",data$SNP[m-1],",Find in ",job)
                                }
                                n_add <- c(n_pos_1,data$POS[m-1],n_wn,data$CHROM[i],foot)
                                if (is.na(data$class[i])){
                                    data$class[i] <- "near"
                                }
                                break
                            }
                        }
                    }
                }
                if (length(data$POS[m-1])>0){
                    if (n_pos_1 !=data$POS[m-1]){
                        
                        near <- rbind(near,n_add)
                        n_add <- c()
                    }
                }
            }
            
        }
    }
    colnames(near) <- c("p1","p2","info","chr","number")
    
    # 删除重复行
    near_new <- data.frame(matrix(ncol =5)) 
    near_new <- near_new[-1,]
    for (i in 1:nrow(near)){
        if (!identical(near$p1[i],near$p2[i])){
            new_add <- c(near$p1[i],near$p2[i],near$info[i],near$chr[i],near$number[i])
            near_new <- rbind(near_new,new_add)
        }
    }
    colnames(near_new) <- c("p1","p2","info","chr","number") 
    near <- near_new
    
    ### 连续标记回帖参考基因组----
    OK <- 0 #成功添加的个数
    for (chr in names(chr_list)){
        
        for (i in 1:nrow(near)){
            if (identical(near$chr[i],chr)){
                my_a <- which.min(abs(as.numeric(near$p1[i]) - as.numeric(chr_list[[chr]]$X3G.Start.1)))
                my_b <- which.min(abs(as.numeric(near$p2[i]) - as.numeric(chr_list[[chr]]$X3G.Start.1)))
                chr_list[[chr]]$out[my_a:my_b] <- near$info[i]
                OK <- OK + (my_b - my_a + 1)
            }
        }
        # cli::cli_alert_success(str_c("[Chromosomes are currently being processed]：",chr))
    }
    num_near <- OK
    cli::cli_alert_success(str_c("[Add success Near Region Number]：",num_near))
    
    ### 迭代添加单标记信息----
    OK <- 0 #成功添加的个数
    for (chr in names(chr_list)){
        
        for (i in 1:nrow(single)){
            if (identical(single$chr[i],chr)){
                index <- which.min(abs(as.numeric(single$positon[i]) - as.numeric(chr_list[[chr]]$X3G.Start.1)))
                chr_list[[chr]]$out[index] <- single$info[i]
                OK <- OK + 1
            }
        }
        # cli::cli_alert_success(str_c("[Chromosomes are currently being processed]：",chr))
    }
    num_single <- OK
    cli::cli_alert_success(str_c("[Add success Single Number]：",num_single))
    ### 迭代添加首尾标记 =====================================================================
    cli::cli_alert_success(str_c("[Run Over]：",Sys.time()))
    
    all_near[[id]] <- near
    all_single[[id]] <- single
    
    write_excel_csv(single,paste0("./out/2_SNP_Infomation/",id,"_single.csv"))
    write_excel_csv(near,paste0("./out/2_SNP_Infomation/",id,"_near.csv"))
    write_excel_csv(data,paste0("./out/2_SNP_Infomation/",id,"_data.csv"),na = "near")
    
    out <- do.call(rbind,lapply(chr_list,function(x)x))
    write_excel_csv(out,paste0("./out/3_Gene_Maping_Result/",job,".snp.csv"),na = "")
    # saveRDS(all_near,file = "./out/SNP_info.rds")
}

# step3：save2write----
id_list_step3 <- list.files("./out/3_Gene_Maping_Result/",pattern = "*.csv")
out_lines <- nrow(out)
tem <- Ref[1:out_lines,1:7] # 需要格外注意：tem和out文件必须一一对应
index <- 8 #从第8列开始标注
# 提取并标注注释信息
for (id in id_list_step3){
    now <- read.csv(paste0("./out/3_Gene_Maping_Result/",id))
    atom <- str_split(id,"[.]")
    job <- atom[[1]][1]
    tem <- bind_cols(tem,now[,8])
    colnames(tem)[index] <- job
    index <- index + 1
    cli::cli_alert_success(str_c("[ Run ]： ",id,Sys.time()))
}

# 将多列信息合并为一列（优化算法）
tem <- tem %>% as.data.frame()
tem$all <- NA
tem_rm_na <- tem[,colSums(is.na(tem)) < nrow(tem)]
tem_rm_na_info <- tem_rm_na[,9:ncol(tem_rm_na)-1]
tem_rm_na$all <- apply(tem_rm_na_info, 1, function(x){
    x <- na.omit(x) # 删除NA值
    x <- x[nchar(x) > 3] # 保留字符长度大于3的元素
    paste(x, collapse = " ; ") # 使用分号作为分隔符连接字符串
 }
)

tem_rm_na <- tem_rm_na[,c(1:7,60,8:59)]

## 结果保存----
tem_rm_na$Num <- apply(tem_rm_na,1,function(x){
    tmp <- x[9:60]
    sum(!is.null(tmp))
})
tem_rm_na$Num <- tem_rm_na$all %>% str_count(";")+1

tem_rm_na <- tem_rm_na[,c(1:8,61,9:60)]

write_tsv(tem_rm_na,str_c("./out/4_Rebind_All_Output/",prefix,"_IT_DS_MLM_CPU.Output.final.tsv"))



