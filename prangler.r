# for internal use

library(shiny)
library(listviewer)
library(DT)
library(RJSONIO)
library(rlist)
library(readxl)
options(shiny.maxRequestSize=100*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Genowis数据结构化校验平台"),
  fluidRow(
    column(fileInput("file","上传50例原始文本表格(必须包含“病理号”列)"),width=3),
    column(h4("步骤：1上传xlsx和json文件；2，搜索病理号，在右侧json区修改结果；3 完毕后点击“text”视图，复制结果至左下文本框；4 点击“开始计算”，看结果。若需要比对两组json，保证右侧上传框为算法版本，左下文本框为医生版本。"),width=6),
    column(fileInput("json","上传UTF-8编码的json文件，且key不为空，即不存在\"\":这种字符串，value不能为null，需要首先线下替换（可以是全量json，会根据病理号自动筛选）",width="160%"),width=3)
    #column(radioButtons("IHC_trans","IHC结构",choices=list('{{}{}}'=F, '[{}{}]'=T), selected=F), width=1)
  ),
  fluidRow(
    column(dataTableOutput("raw"),width = 9),
    column(jsoneditOutput("json_result",height = "650px"),width = 3)),
  fluidRow(
    column(textAreaInput("json_final","输入最终的json文件,其中免疫组化为[{}{}]",width='700px'),width = 5),
    column(radioButtons("partial","部分数据",choices=list(全量=F, 状态码1=T),selected=T),width=1),
    column(h5("重复病理号（存在重复病理号无法进行比对）"),textOutput("id_same"),width = 3),
    column(tableOutput("results"),h5("若完整性不为1但漏提病例数为0，是标准版本指标变少（较罕见）造成"),width=3)),
  
  h4("校验细节"),
  column(width=12,
         dataTableOutput("detail"))
  
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  rawFile=reactive({
    rawPath=input$file$datapath
    read_xlsx(rawPath)
  })
  
  transformIHC=function(x){
    part=names(x$emr_info) %>% .[.!="诊断结论"]
    for (i in part){
      check=x$emr_info[[i]]
      if (!is.element("免疫组化",names(check))) {
        next
      }
      else {
        cand=check$免疫组化
        names(cand)=0:(length(cand)-1)
        x$emr_info[[i]]$免疫组化=cand
      }
    }
    return(x)
  }
  
  jsonFile=reactive({
    jsonPath=input$json$datapath
    jsonData=RJSONIO::fromJSON(jsonPath)
    if (is.null(input$file$datapath)){
      raw=jsonData
    }
    else {
      PID=rawFile()$病理号
      raw=list.filter(jsonData,base_info[["病理号"]] %in% PID)
    }	
    lapply(raw, transformIHC)
  })
  
  output$raw=renderDataTable({
    datatable(rawFile(), extensions = list('Buttons',"Scroller"),
              options = list(pageLength = 50,dom = 'Bfrtp',buttons = I('colvis'),scrollY = 550,scroller = TRUE))
  })
  output$json_result=renderJsonedit({
    jsonedit(jsonFile(),height=600,modes=c("text","tree"))
  })
  
  output$id_same=renderText({
    PID=unlist(list.select(jsonFile(),base_info[['病理号']])) 
    paste0(unique(PID[duplicated(PID)]),collapse="@")
  }) 
  
  json_Inspect=function(x,y){
    # x,y are the json files with same length (observation lines), y is the final result, or standard set
    PID=unlist(list.select(y,base_info[['病理号']]))
    #PID=rawFile()$病理号
    result=NULL
    recall_count=list(orig=0,missing=0)
    precision_count=list(orig=0,incorrect=0)
    withProgress(message = 'JSON校验中', value = 0, {
      for (id in PID) {
        final=list.filter(y, base_info[["病理号"]]==id) %>% unlist()%>% t()%>% as.data.frame(stringsAsFactors=F)
        original=list.filter(x, base_info[["病理号"]]==id) 
        if (length(original)==0){
          recall_count=list(orig=recall_count$orig+ncol(final),missing=recall_count$missing+ncol(final))
          recall_eval=paste0(colnames(final),collapse = "@")
          precision_count=precision_count
          precision_eval=precision_eval
        }
        else {
          original= original %>% unlist()%>% t()%>% as.data.frame(stringsAsFactors=F)
          missing_keys=setdiff(colnames(final),colnames(original)) 
          recall_count=list(orig=recall_count$orig+ncol(final),missing=recall_count$missing+length(missing_keys))
          recall_eval=paste0(missing_keys,collapse = "@")
          common_keys=intersect(colnames(final),colnames(original)) 
          final_common=final[common_keys]
          original_common=original[common_keys]
          diff_values_idx=which(t(final_common)!=t(original_common))  
          final_values=final_common[diff_values_idx]
          original_values=original_common[diff_values_idx]
          precision_count=list(orig=precision_count$orig+length(common_keys),incorrect=precision_count$incorrect+length(diff_values_idx))	  
          precision_eval=paste(common_keys[diff_values_idx],"@医生：",final_values, "@算法：",original_values,sep="",collapse = "♥")	      
        }
        result=data.frame(PID=id,完整性评估 =recall_eval,准确性评估=precision_eval,stringsAsFactors = F) %>% rbind(result,.) %>% subset(.,nchar(准确性评估)+nchar(完整性评估)>8)
        #setTxtProgressBar(txtProgressBar(min=0,max=1,style=3),value=which(PID==id)/length(PID))
        incProgress(1/length(PID),detail = paste("正在计算第", which(PID==id),"个",sep=""))
      }
    })    
    recall_overall=1-(recall_count$missing/recall_count$orig)
    precision_overall=1-(precision_count$incorrect/precision_count$orig)
    return(list(完整性=recall_overall, 准确性=precision_overall,结果=result))
  }
  
  result_all=reactive({
    ready=fromJSON(input$json_final) %>% lapply(., transformIHC)
    if (input$partial==F){
      ready=ready}
    else {	
      ready= list.filter(ready,base_info[["状态码"]]==1)
    }
    json_Inspect(jsonFile(),ready)
  })  
  output$results=renderTable({
    data.frame(完整性=result_all()$完整性,准确性=result_all()$准确性,出错病例数=sum(nchar(result_all()$结果$准确性评估)>8),漏提病例数=sum(nchar(result_all()$结果$完整性评估)>0))	
  },digits=6)
  
  output$detail=renderDataTable({
    datatable(result_all()$结果, rownames=F,filter = 'top',options = list(pageLength = 25))    
  })
}

# Run the application 
shinyApp(ui = ui, server = server,options=list(port = 0331,host = "0.0.0.0",launch.browser=F, quiet=T))

