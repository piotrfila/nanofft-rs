#![no_std]
#![no_main]

extern crate panic_usb_boot;

#[rtic::app(
    device = rp_pico::hal::pac,
    dispatchers = [TIMER_IRQ_1]
)]
mod app {
    use rp2040_monotonic::Rp2040Monotonic;

    use rp_pico::{
        hal,
        hal::gpio::bank0::Gpio25,
        Pins,
        XOSC_CRYSTAL_FREQ,
    };

    use core::mem::MaybeUninit;
    use embedded_hal::digital::v2::{OutputPin, ToggleableOutputPin};
    // use panic_usb_boot as _;
    use usb_device::{class_prelude::*, prelude::*};

    const MONO_NUM: u32 = 1;
    const MONO_DENOM: u32 = 1_000_000;
    type Duration = fugit::Duration::<u64, MONO_NUM, MONO_DENOM>;
    const LED_TOGGLE_DELAY: Duration = Duration::from_ticks(200_000); // 200ms

    #[monotonic(binds = TIMER_IRQ_0, default = true)]
    type Rp2040Mono = Rp2040Monotonic;

    #[shared]
    struct Shared { }

    #[local]
    struct Local {
        fifo: rp_pico::hal::sio::SioFifo,
        led: hal::gpio::Pin<Gpio25, hal::gpio::FunctionSio<hal::gpio::SioOutput>, hal::gpio::PullDown>,
        serial: usbd_serial::SerialPort<'static, hal::usb::UsbBus>,
        usb_dev: UsbDevice<'static, hal::usb::UsbBus>
    }

    #[init(local=[
        // Task local initialized resources are static
        // Here we use MaybeUninit to allow for initialization in init()
        // This enables its usage in driver initialization
        usb_bus: MaybeUninit<UsbBusAllocator<hal::usb::UsbBus>> = MaybeUninit::uninit(),
    ])]
    fn init(ctx: init::Context) -> (Shared, Local, init::Monotonics) {
        let mut pac = ctx.device;

        // Configure the clocks, watchdog - The default is to generate a 125 MHz system clock
        let mut watchdog = hal::watchdog::Watchdog::new(pac.WATCHDOG);
        let clocks = hal::clocks::init_clocks_and_plls(
            XOSC_CRYSTAL_FREQ,
            pac.XOSC,
            pac.CLOCKS,
            pac.PLL_SYS,
            pac.PLL_USB,
            &mut pac.RESETS,
            &mut watchdog,
        )
        .ok()
        .unwrap();

        // Init LED pin
        let mut sio = hal::Sio::new(pac.SIO);
        let gpio = Pins::new(
            pac.IO_BANK0,
            pac.PADS_BANK0,
            sio.gpio_bank0,
            &mut pac.RESETS,
        );
        let mut led = gpio.led.into_push_pull_output();
        led.set_high().unwrap();

        let mono = Rp2040Mono::new(pac.TIMER);

        heartbeat::spawn_after(LED_TOGGLE_DELAY).unwrap();

        let mut mc = hal::multicore::Multicore::new(
            &mut pac.PSM,
            &mut pac.PPB,
            &mut sio.fifo,
        );
        let core1 = &mut mc.cores()[1];
        let core1_stack: &mut MaybeUninit<hal::multicore::Stack<1024>> = unsafe { &mut *(0x20040000 as *mut _) };
        let core1_stack = core1_stack.write(hal::multicore::Stack::new());
        let _ = core1.spawn(&mut core1_stack.mem, move || {
            let pac = unsafe { hal::pac::Peripherals::steal() };
            let mut sio = hal::Sio::new(pac.SIO);

            let buf: &mut [(i16, i16); 65536] = unsafe { &mut *(0x20000000 as *mut _) };

            let max_size = 16384;

            loop {
                let mut size: usize = 4;

                let _ = sio.fifo.read_blocking();

                while size <= max_size {
                    let repeats = max_size as u32 * max_size.trailing_zeros() / (size as u32 * size.trailing_zeros());
                    // let repeats = max_size / size;
                    let buf = &mut buf[..size];

                    let t1 = crate::app::monotonics::now();
                    
                    for _ in 0..repeats {
                        let _ = nanofft::i16::fft_pairs_dyn(x);
                    }
                    let t2 = crate::app::monotonics::now();

                    let time_in_millionths_of_cycle = (t2.ticks().wrapping_sub(t1.ticks()) as u64) * 125_000_000;
                    let time_per_iteration = time_in_millionths_of_cycle / (repeats as u64);
                    let time_per_multiply = time_per_iteration / ((size as u32 * size.trailing_zeros()) as u64);
                    
                    sio.fifo.write_blocking(size as _);
                    sio.fifo.write_blocking(time_per_multiply.try_into().unwrap());
                    size <<= 1;
                }
                sio.fifo.write_blocking(0);
                sio.fifo.write_blocking(0);

            }
        });

        let fifo = sio.fifo;

        let usb_bus = ctx.local.usb_bus.write(
            UsbBusAllocator::new(hal::usb::UsbBus::new(
                pac.USBCTRL_REGS,
                pac.USBCTRL_DPRAM,
                clocks.usb_clock,
                true,
                &mut pac.RESETS,
            ))
        );

        let serial = usbd_serial::SerialPort::new(usb_bus);

        let usb_dev = UsbDeviceBuilder::new(usb_bus, UsbVidPid(0x16c0, 0x27dd))
            .manufacturer("Fake company")
            .product("Serial port")
            .serial_number("TEST")
            .device_class(2) // from: https://www.usb.org/defined-class-codes
            .build();

        // Return resources and timer
        (
            Shared {},
            Local { fifo, led, serial, usb_dev },
            init::Monotonics(mono),
        )
    }

    #[task(local = [led])]
    fn heartbeat(ctx: heartbeat::Context) {
        _ = ctx.local.led.toggle();
        heartbeat::spawn_after(LED_TOGGLE_DELAY).unwrap();
    }

    #[task(binds = USBCTRL_IRQ, local = [fifo, serial, usb_dev])]
    fn usb_poll(ctx: usb_poll::Context) {
        let fifo = ctx.local.fifo;
        let serial = ctx.local.serial;
        let usb_dev = ctx.local.usb_dev;

        if usb_dev.poll(&mut [serial]) {
            let mut buf = [0u8; 16];
            match serial.read(&mut buf) {
                Ok(0) | Err(_) => { }
                Ok(count) => {
                    for i in 0..count {
                        if buf[i] == 'r' as u8 {
                            hal::rom_data::reset_to_usb_boot(0, 0);
                        }
                    }
                    let _ = fifo.write(0);
                }
            }
        }

        if let Some(mut ret) = fifo.read() {
            let mut buf = [' ' as u8; 16];
            buf[0] = '.' as u8;
            buf[buf.len() - 2] = '\r' as u8;
            buf[buf.len() - 1] = '\n' as u8;
            let mut wr_ptr = &mut buf[..];
            
            for i in (0..(wr_ptr.len() - 2)).rev() {
                if i % 4 == 2 { continue }
                wr_ptr[i] = '0' as u8 + ((ret % 10) as u8);
                ret /= 10;
                if ret == 0 {
                    wr_ptr = &mut wr_ptr[i..];
                    break;
                }
            }

            let mut wr_ptr = &buf[..];

            while !wr_ptr.is_empty() {
                let _ = serial.write(wr_ptr).map(|len| {
                    wr_ptr = &wr_ptr[len..];
                });
            }
        }
    }
}