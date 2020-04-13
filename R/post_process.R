post_process=function(se,percent=0.9){

  #post process the enhanced network
  diag(se)=0
  q=quantile(se,percent)

  #denoise the network
  se[se<q]=0
  se[se>q]=1

  #make the network symmetric
  se=(se+t(se))/2
  se[se==0.5]=1

  return(se)
}
