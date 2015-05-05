  cl_int status;
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  std::vector<cl::Device> devices;
  platforms[1].getDevices(CL_DEVICE_TYPE_ALL, &devices);

  cl::Context context(devices);
  cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);
  cl::Program program = getProgram(context, "sw_f32.cl");
  status = program.build(devices);
  std::string msg;
  program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &msg);
  std::cout << msg << std::endl;
  std::cout << "program built" << std::endl;
  
  cl::Buffer buffer_sm(context, CL_MEM_READ_ONLY, (128 * 128) * sizeof(float));
  queue.enqueueWriteBuffer(buffer_sm, CL_TRUE, 0, (128 * 128) * sizeof(float), sm);

  cl::NDRange global(buffer_width);
  //cl::NDRange local(1024);
  
  
  
  /*
  std::vector<cl::Event> write_buffer_events(buffer_count * 2);
  cl::Kernel sw_f32(program, "sw_f32");
  sw_f32.setArg(2, buffer_sm);
  sw_f32.setArg(3, gap_open);
  sw_f32.setArg(4, gap_extend);
*/
  
  
  
  while(batches[i]->fill(reader1, reader2) > 0)
  {
    b++;
    /*
    std::vector<cl::Event> write_buffer_events(2);
    queue.enqueueWriteBuffer(batches[i]->buffer_1, CL_FALSE, 0,
                             buffer_width * batches[i]->max_length1, batches[i]->seqs1.data(),
                             NULL, &write_buffer_events[0]);
    queue.enqueueWriteBuffer(batches[i]->buffer_2, CL_FALSE, 0,
                             buffer_width * batches[i]->max_length2, batches[i]->seqs2.data(),
                             NULL, &write_buffer_events[1]);
    sw_f32.setArg(0,  batches[i]->buffer_1);
    sw_f32.setArg(1,  batches[i]->buffer_2);
    sw_f32.setArg(5,  batches[i]->buffer_flags);
    sw_f32.setArg(6,  batches[i]->buffer_max_scores);
    sw_f32.setArg(7,  batches[i]->buffer_i_max);
    sw_f32.setArg(8,  batches[i]->buffer_j_max);
    sw_f32.setArg(9,  batches[i]->max_length1 + 1);
    sw_f32.setArg(10, batches[i]->max_length2 + 1);

    cl::Event kernel_event;
    status = queue.enqueueNDRangeKernel(sw_f32, cl::NullRange, global, cl::NullRange,
    &write_buffer_events, &kernel_event);

    if(status != CL_SUCCESS)
    {
      std::cout << "Error " << status << std::endl;
    }

    std::vector<cl::Event> kernel_events;
    kernel_events.push_back(kernel_event);
    
    queue.enqueueReadBuffer(batches[i]->buffer_max_scores, CL_FALSE, 0,
    buffer_width * sizeof(float), batches[i]->max_scores.data(), &kernel_events);
    //queue.enqueueReadBuffer(batches[i]->buffer_flags, CL_FALSE, 0,
    //buffer_length * buffer_length * buffer_width, batches[i]->flags, &kernel_events);
    
    for(int j =0;j<1;j++)
      std::cout <<  batches[i].max_scores[j]<< std::endl;
    */
    i++;
    i %= buffer_count;
    std::cout << "Batch "<< b << " processed" << std::endl;    
  }
  
  
    for(int i = 0; i < buffer_count; i++)
  {
#ifdef OCL_SUPPORT
    batches.push_back(new aln::AlnBatch(buffer_width, buffer_length, context));
#else
    batches.push_back(new aln::AlnBatch(buffer_width, buffer_length));
#endif
  }
